import math
import os
import re
import sys
import shutil
import asyncio
import pathlib
import logging
import tempfile
import argparse
import itertools
import subprocess

import gmpy2
from aiofiles import os as aios


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def locate_gpu_ecm_install():
    install_location = shutil.which("ecm")
    result = subprocess.run([install_location, "-printconfig"], capture_output=True)
    compile_options = result.stdout.decode("utf-8")
    if "WITH_GPU = 1" in compile_options:
        return install_location
    eprint("Could not find GPU ECM binary location, maybe try `make install` in your ecm source dir")
    sys.exit(1)


def determine_optimal_num_gpu_curves():
    return 8192

def get_b1_b2_plan():
    # replace this with dynamic b1 level choice
    return [
        (9562627, 531354383),
        (19673032, 1445101698),
        (19672032, 1444962986),
        (24852366, 2163537440),
        (46502049, 5757597518),
        (46501349, 5757473846),
        (46735326, 5798811809),
        (91219077, 17211831452),
        (109526746, 23522817677),
        (110359262, 23807037137),
        (110356563, 23806115703),
        (211471250, 70463572846),
        (265889508, 97002133243),
        (265866810, 96991063939),
        (265053791, 96594572860),
        (265065390, 96600229432),
        (265049392, 96592427567),
        (356693127, 144620131535),
        (356702826, 144625247257),
        (356705225, 144626512606),
    ]


# async task to send all the inputs, flush, and then close the stdin, i feel like this shouldn't be so difficult
async def handle_stdin(lines, stdin):
    # try:
    stdin.write("\n".join(lines).encode())
    await stdin.drain()
    stdin.close()
    # except BaseException as e:
    #     eprint(e)


async def kill_processes(procs):
    for proc in procs:
        try:
            if proc is not None:
                await proc.kill()
        except (OSError, TypeError) as e:
            logging.debug(e)

def cancel_tasks(tasks):
    for task in tasks:
        task.cancel()


# transform list of factors into all coprime terms
# https://cr.yp.to/lineartime/dcba-20040404.pdf
# or, find coprime base of the factors
def reduce_to_coprimes(factors):
    #TODO
    pass

# new factors were found, calculate remaining cofactor and if we are done
def new_factors_found(input_number, old_cofactor, found_factors):
    # found factors sometimes includes a composite that includes one of the other factors, due to gpu stage 1

    # test them all for primality (should be quick with ecm found factor size)
    prime_found_factors = set()
    composite_found_factors = set()
    for found_factor in found_factors:
        if gmpy2.is_prime(found_factor):
            prime_found_factors.add(found_factor)
        else:
            composite_found_factors.add(found_factor)

    # check if any composite factors contain new prime factors not in prime found factors
    any_found = True
    while any_found:
        any_found = False
        for composite_found_factor in composite_found_factors:
            for prime_found_factor in prime_found_factors:
                if gmpy2.is_divisible(composite_found_factor, prime_found_factor):
                    new_cofactor = composite_found_factor // prime_found_factor
                    any_found = True
                    composite_found_factors.remove(composite_found_factor)
                    if gmpy2.is_prime(new_cofactor):
                        if new_cofactor not in prime_found_factors:
                            prime_found_factors.add(new_cofactor)
                    else:
                        if new_cofactor not in composite_found_factors:
                            composite_found_factors.add(new_cofactor)
                    break
            if any_found:
                break
    found_factors = list(prime_found_factors) + list(composite_found_factors)

    # if all found_factors recreate the input_number, we're done
    # it's possible there are composites in found_factors, but if ECM found them,
    # it's fine to leave them in the factors I guess
    factor_prod = math.prod(found_factors)
    if input_number == factor_prod:
        return 1, found_factors, True
    cofactor, remainder = gmpy2.f_divmod(input_number, factor_prod)
    assert remainder == 0, f"Expected {input_number} % {factor_prod} == 0\n{found_factors}"
    # could return here, but check for prime powers in the factors in case we missed one
    for found_factor in set(found_factors):
        while gmpy2.is_divisible(cofactor, found_factor):
            # add extra prime factor
            found_factors.append(found_factor)
            cofactor = cofactor // found_factor

    if gmpy2.is_prime(cofactor):
        found_factors.append(int(cofactor))
        cofactor = 1

    return cofactor, found_factors, cofactor == 1





async def main():
    input_numbers = []
    num_threads = 16
    curves_per_batch = determine_optimal_num_gpu_curves()
    param = 3
    exit_after_one = True
    with tempfile.TemporaryDirectory() as tmpdir:
        if not sys.stdin.isatty():
            input_numbers = list(map(int, sys.stdin.read().strip().split("\n")))
        install_location = locate_gpu_ecm_install()
        gpu_proc = None
        old_gpu_proc = None
        for input_number in input_numbers:
            if gpu_proc is not None:
                await kill_processes([gpu_proc])
            if old_gpu_proc is not None:
                await kill_processes([old_gpu_proc])
            gpu_proc = None
            old_gpu_proc = None
            eprint(f"# N = {input_number}")
            fully_factored = gmpy2.is_prime(input_number)
            cofactor = input_number
            plan = get_b1_b2_plan()
            temp_save_file_path = None
            found_factors = []
            for i in range(len(plan)):
                if fully_factored or (found_factors and exit_after_one):
                    print(f"{input_number}={found_factors}")
                    break
                b1, b2 = plan[i]
                old_b1, old_b2 = plan[i-1] if i > 0 else plan[i]
                if gpu_proc:  # wait for old gpu_proc to finish before continuing with next
                    # print 0 curves as status to show GPU ecm is running
                    # eprint(f" waiting for gpu_proc", end="\r")
                    eprint(f" {0: >{len(str(curves_per_batch))}}/{curves_per_batch}@{old_b1},{old_b2},{param}  ", end="\r")
                    # TODO add progress bar:
                    # Computing 11366 bits/call, 681174/10035736 (6.8%), ETA 73 + 5 = 78 seconds (~10 ms/curves)
                    # Computing 11366 bits/call, 1817774/10035736 (18.1%), ETA 63 + 14 = 78 seconds (~9 ms/curves)
                    # Computing 11366 bits/call, 2954374/10035736 (29.4%), ETA 55 + 23 = 77 seconds (~9 ms/curves)
                    # Computing 11366 bits/call, 4090974/10035736 (40.8%), ETA 46 + 32 = 77 seconds (~9 ms/curves)
                    outs, errs = await gpu_proc.communicate()
                    gpu_return = gpu_proc.returncode
                    outs = outs.decode() if outs else ""
                    errs = errs.decode() if errs else ""
                    sigma_match = re.search(r"sigma=(\d+):(\d+)-(?:\d+:)?(\d+)", outs)
                    assert sigma_match is not None, "Couldn't find sigma range in ecm gpu output"
                    param = int(sigma_match.group(1))
                    low_sigma = int(sigma_match.group(2))
                    high_sigma = int(sigma_match.group(3))

                    factor_matches = re.findall(
                        r"GPU: factor (\d+) found in Step 1 with curve (\d+) \(-sigma (\d+):(\d+)\)", outs)
                    # save earliest hit per factor, (factor, curve_n, param, sigma)
                    gpu_factors = dict(map(lambda x: (int(x[0]), tuple(map(int, x[1:]))), reversed(factor_matches)))
                    if gpu_factors:
                        found_factors.extend(gpu_factors.keys())
                        cofactor, found_factors, fully_factored = new_factors_found(input_number, cofactor, found_factors)
                    if fully_factored or (found_factors and exit_after_one):
                        eprint()
                        print(f"{input_number}={found_factors}")
                        break

                    # eprint(outs)
                    # eprint(errs)
                    # print(f"\n\nsigma={param}:{low_sigma}-{param}:{high_sigma}")
                    if gpu_return == 8:  # found factor(s)
                        eprint("Found factors, but couldn't find in stdout, weird...")
                        # eprint(outs)
                        # eprint(errs)
                        sys.exit(1)
                    elif gpu_return == 0:  # no factors found, but no errors
                        pass
                    else:  # errors probably
                        eprint(errs)
                        raise Exception(f"{errs}\nGPU process exited with {gpu_return}")
                old_temp_save_file_path = temp_save_file_path
                old_gpu_proc = gpu_proc
                temp_save_file_name = f"{hash(input_number)}_{b1}_stg1_residues"
                temp_save_file_path = os.path.join(tmpdir, temp_save_file_name)
                gpu_ecm_args = [install_location, "-gpu", "-v", "-save", temp_save_file_path, str(b1), str(0)]
                gpu_proc = await asyncio.subprocess.create_subprocess_exec(*gpu_ecm_args,
                                                                           stdout=asyncio.subprocess.PIPE,
                                                                           stdin=asyncio.subprocess.PIPE,
                                                                           stderr=asyncio.subprocess.PIPE)
                # synchronously write, flush, and close the one number to gpu ecm, shouldn't block
                gpu_proc.stdin.write(str(cofactor).encode())
                await gpu_proc.stdin.drain()
                gpu_proc.stdin.close()


                # get old gpu return val, or skip this loop if it's first time and it doesn't exist
                # old gpu proc should already be finished by the time we make the new gpu proc
                if old_gpu_proc:
                    gpu_return = await old_gpu_proc.wait()
                else:
                    continue

                if gpu_return == 8:
                    # factors found return code, should just go to other loop and get them
                    # running stage-2 on residues likely just reveals same factors
                    eprint("\nstage-1 factor found\n")
                    continue

                if old_temp_save_file_path and await aios.path.isfile(old_temp_save_file_path):
                    # run cpu curves
                    residue_lines = pathlib.Path(old_temp_save_file_path).read_text().strip().split("\n")
                    enumerated_residue_lines = sorted(map(lambda x: (x[0] % num_threads, x[1]),
                                               enumerate(
                                                   residue_lines)),
                                           key=lambda x: x[0])
                    grouped_residues = itertools.groupby(enumerated_residue_lines, key=lambda x: x[0] % num_threads)
                    cpu_ecm_args = [install_location, "-resume", "-", str(old_b1), str(old_b2)]
                    cpu_procs = []
                    cpu_proc_residues = {}
                    tasks = []
                    for i, residues in grouped_residues:
                        cpu_proc_residues[i] = list(map(lambda x: x[1], residues))
                        cpu_proc = await asyncio.create_subprocess_exec(*cpu_ecm_args,
                                                                        stdout=asyncio.subprocess.PIPE,
                                                                        stdin=asyncio.subprocess.PIPE,
                                                                        stderr=asyncio.subprocess.STDOUT)
                        tasks.append(asyncio.get_event_loop().create_task(handle_stdin(cpu_proc_residues[i], cpu_proc.stdin)))
                        cpu_procs.append(cpu_proc)
                    all_cpu_processes_done = False
                    count = 0
                    total_curves = len(residue_lines)
                    total_curves_digits = len(str(total_curves))
                    while not all_cpu_processes_done:
                        if found_factors:
                            # kill cpu processes now that we have factor
                            cancel_tasks(tasks)
                            await kill_processes(cpu_procs)
                            cofactor, found_factors, fully_factored = new_factors_found(input_number, cofactor, found_factors)
                            # break out of CPU loop to go see if GPU is done
                            break

                        all_cpu_processes_done = True
                        for cpu_proc in cpu_procs:
                            line = await cpu_proc.stdout.readline()
                            if line:
                                line = line.decode().strip()
                                all_cpu_processes_done = False
                                if not line.startswith("Step 1") and not line.startswith("Using") and not line.startswith("Input number is") and not line.startswith("Resuming ECM") and not line.startswith("GMP-ECM"):
                                    if line.startswith("Step 2"):
                                        count += 1
                                        eprint(f" {count: >{total_curves_digits}}/{total_curves}@{old_b1},{old_b2},{param}  ", end="\r")
                                    elif line.startswith("********** Factor found in step 2:"):
                                        found_factors.append(int(line.strip().split(" ")[-1]))
                                        break
                                    else:
                                        eprint(f"\n{line}\n")
                                        eprint(f"sys exit")
                                        sys.exit(1)  # not sure how to parse yet

                    eprint()
        await asyncio.sleep(1)
        if gpu_proc is not None:
            await kill_processes([gpu_proc])
        if old_gpu_proc is not None:
            await kill_processes([old_gpu_proc])
        # this is necessary for some godforsaken reason, otherwise there's an asyncio exception from garbage collection
        await asyncio.sleep(0.1)



if __name__ == "__main__":
    asyncio.run(main(), debug=False)
