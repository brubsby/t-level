import math
import os
import re
import sys
import shutil
import signal
import asyncio
import pathlib
import logging
import tempfile
import argparse
import itertools
import subprocess

import gmpy2
from aiofiles import os as aios

interrupt_level = 0
gpu_proc = old_gpu_proc = cpu_procs = tasks = None
loop = asyncio.new_event_loop()
message = ""

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def positive_integer(arg):
    val = int(arg)
    if val < 1:
        raise ValueError(f"{arg} not a positive integer")
    return val


def nonnegative_integer(arg):
    val = int(arg)
    if val < 0:
        raise ValueError(f"{arg} not a non-negative integer")
    return val


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
        (100, 600000),
        (136808, 8000000),
        (500000, 15000000),
        (901309, 21288716),
        (1669233, 52732403),
        (2994200, 131737403),
        (5026497, 352580564),
        (6872112, 531505103),
        (14019466, 1443854888),
        (17829816, 2185199673),
        (27253419, 4313913596),
        (33115351, 5777998123),
        (33117865, 5778623383),
        (77871012, 23621754325),
        (77606595, 23492099779),
        (77744024, 23559486878),
        (149939404, 70874529876),
        (150583940, 71319614527),
        (187054292, 96504234003),
        (186925087, 96415011427),
        (186928084, 96417081007),
        (251865710, 144528548178),
        (251850112, 144516923636),
        (319326164, 194825966526),
        (316463050, 192688956109),
        (524397977, 388625202894),
    ]

# def get_b1_b2_plan():
#     # replace this with dynamic b1 level choice
#     return [
#         (100, 100000000),
#         (40434702, 23524201837),
#         (59868512, 46897505668),
#         (77569142, 70713526451),
#         (96980101, 96830777812),
#         (97044194, 96917014303),
#         (130465452, 144878169431),
#         (163525145, 192530863837),
#         (163516546, 192518469121),
#         (271217975, 388385545690),
#         (271227974, 388412217655),
#         (272381659, 391489629978),
#         (581431351, 1590060278137),
#         (690639929, 2384784962810),
#         (690623430, 2384664897493),
#         (582285965, 1596279412697),
#         (821583233, 3218422750288),
#         (1172685119, 6472256728113),
#         (1442021083, 9709869363672),
#         (1163349953, 6371594488808),
#         (1163357152, 6371671817551),
#         (1428911694, 9537400540544),
#     ]


# async task to send all the inputs, flush, and then close the stdin, i feel like this shouldn't be so difficult
async def handle_stdin(queue, stdin):
    try:
        while not queue.empty():
            to_send = (await queue.get() + "\n")
            # eprint(f"sending {to_send}", end="")
            stdin.write(to_send.encode())
            await stdin.drain()
            queue.task_done()
        stdin.close()
    # handle case where process dies or is canceled
    except (BrokenPipeError, ConnectionResetError):
        pass
    except Exception as e:
        # eprint()
        # eprint(e)
        raise e


async def kill_processes(procs, log=True):
    for proc in procs:
        try:
            if proc is not None:
                await proc.kill()
        except (OSError, TypeError) as e:
            if log:
                logging.debug(e)


async def kill_gpu_procs(log=True):
    global gpu_proc, old_gpu_proc
    if gpu_proc is not None:
        await kill_processes([gpu_proc], log)
    if old_gpu_proc is not None:
        await kill_processes([old_gpu_proc], log)
    gpu_proc = None


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


async def shutdown(errcode):
    eprint()
    global gpu_proc, old_gpu_proc, cpu_procs, tasks
    if gpu_proc is not None:
        await kill_processes([gpu_proc], log=False)
    if old_gpu_proc is not None:
        await kill_processes([old_gpu_proc], log=False)
    if cpu_procs is not None:
        await kill_processes(cpu_procs, log=False)
    if tasks is not None:
        for task in tasks:
            try:
                await task
            except BrokenPipeError:
                pass

    await asyncio.sleep(0.1)
    loop.stop()


async def main():
    global cpu_procs, gpu_proc, old_gpu_proc, tasks, interrupt_level, message

    __version__ = "0.0.2"

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"")
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="verbosity (-v, -vv, etc)")
    # parser.add_argument(
    #     "-q",
    #     type=str,
    #     action="store",
    #     dest="expression",
    #     help="direct expression input (currently only decimal expansions supported)",
    # )
    # parser.add_argument(
    #     "-i",
    #     "--input",
    #     type=str,
    #     action="store",
    #     dest="filename",
    #     help="file containing composites to run ecm on",
    # )
    # parser.add_argument(
    #     "-w",
    #     "--work",
    #     action="store",
    #     dest="work",
    #     type=float,
    #     help="existing t-level of work done, determines starting point of work"
    # )
    # parser.add_argument(
    #     "-p",
    #     "--pretest",
    #     action="store",
    #     dest="pretest",
    #     type=float,
    #     help="the desired t-level to reach in deceimal digits, quits after reaching"
    # )
    # parser.add_argument(
    #     "-t",
    #     "--threads",
    #     action="store",
    #     dest="threads",
    #     type=positive_integer,
    #     default=os.cpu_count(),
    #     help="number of threads to use in CPU stage-2"
    # )
    # parser.add_argument(
    #     "--gpudevice",
    #     action="store",
    #     dest="gpu_device",
    #     type=nonnegative_integer,
    #     help="use device <GPU_DEVICE> to execute GPU code (by default, CUDA chooses)"
    # )
    # parser.add_argument(
    #     "--gpucurves",
    #     action="store",
    #     dest="gpu_curves",
    #     type=positive_integer,
    #     help="compute on <GPU_CURVES> curves in parallel on the GPU (by default, CUDA chooses)"
    # )
    # parser.add_argument(
    #     "-o",
    #     "--one",
    #     action="store",
    #     dest="exit_after_one",
    #     help="stop ECM curves on composite after one factor is found, equivalent to -x 1"
    # )
    # parser.add_argument(
    #     "-x",
    #     "--exitafter",
    #     action="store_true",
    #     dest="exit_after",
    #     type=positive_integer,
    #     help="stop ECM curves on a composite after <EXIT_AFTER> factors are found, or cofactor is prime"
    # )
    # parser.add_argument(
    #     "--tune",
    #     action="store_true",
    #     dest="tune",
    #     help="run tuning"
    # )

    args = parser.parse_args()

    loglevel = logging.WARNING
    if args.verbose > 0:
        loglevel = logging.INFO
    if args.verbose > 1:
        loglevel = logging.DEBUG
    logging.basicConfig(level=loglevel, format="%(message)s")

    async def handle_signals(signame):
        global interrupt_level, message
        if signame == 'SIGINT':
             interrupt_level += 1
        if signame == 'SIGTERM':
            interrupt_level += 10
        inner_message = ""
        if interrupt_level == 1:
            inner_message = "won't start another b-level with GPU..."
        elif interrupt_level == 2:
            inner_message = "killing gpu process, finishing CPU curves..."
            await kill_gpu_procs()
        elif interrupt_level == 3:
            inner_message = "killing all remaining processes and quitting..."
        message = f"# Interrupt level: {interrupt_level}{', ' + inner_message if inner_message else ''}"
        # eprint(message, end="\r")
        if interrupt_level >= 4:
            eprint("\n# Shutting Down...", end="")
            await shutdown(1)

    for signame in ('SIGINT', 'SIGTERM'):
        loop.add_signal_handler(getattr(signal, signame),
                                lambda signame=signame: asyncio.create_task(handle_signals(signame)))

    input_numbers = []
    factor_size_sigma_report_threshold = 60
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
            await kill_gpu_procs(True)
            gpu_proc = None
            old_gpu_proc = None
            eprint(f"# N = {input_number}", end="")
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
                    eprint()
                    # print 0 curves as status to show GPU ecm is running
                    # eprint(f" waiting for gpu_proc", end="\r")
                    eprint(f" {0: >{len(str(curves_per_batch))}}/{curves_per_batch}@{old_b1},{old_b2},{param}  {message}", end="\r")
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
                    try:
                        assert sigma_match is not None, "Couldn't find sigma range in ecm gpu output"
                        param = int(sigma_match.group(1))
                        low_sigma = int(sigma_match.group(2))
                        high_sigma = int(sigma_match.group(3))
                    except (AssertionError, AttributeError) as e:
                        if not interrupt_level:
                            raise e

                    factor_matches = re.findall(
                        r"GPU: factor (\d+) found in Step 1 with curve (\d+) \(-sigma (\d+):(\d+)\)", outs)
                    # save earliest hit per factor, (factor, curve_n, param, sigma)
                    gpu_factors = dict(map(lambda x: (int(x[0]), tuple(map(int, x[1:]))), reversed(factor_matches)))
                    if gpu_factors:
                        found_factors.extend(gpu_factors.keys())
                        for factor, params in gpu_factors.items():
                            if math.log10(factor) >= factor_size_sigma_report_threshold:
                                eprint(f"********** BIG GPU ECM STAGE-1 HIT: TELL YOUR FRIENDS! SIGMA={params[1]}:{params[2]} **********")
                        cofactor, found_factors, fully_factored = new_factors_found(input_number, cofactor, found_factors)
                        eprint(factor_matches)
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
                        await shutdown(1)
                    elif gpu_return in [0, 143]:  # no factors found, but no errors, 143 is sigterm
                        pass
                    elif interrupt_level >= 2:  # user requested killing GPU, it won't have any results
                        pass
                    else:  # errors probably
                        eprint(errs)
                        eprint(outs)
                        raise Exception(f"{errs}\nGPU process exited with {gpu_return}")
                old_temp_save_file_path = temp_save_file_path
                old_gpu_proc = gpu_proc
                temp_save_file_path = None
                gpu_proc = None
                if interrupt_level == 0:  # only make gpu process if we're not interrupted
                    temp_save_file_name = f"{hash(input_number)}_{b1}_stg1_residues"
                    temp_save_file_path = os.path.join(tmpdir, temp_save_file_name)
                    gpu_ecm_args = [install_location, "-gpu", "-v", "-save", temp_save_file_path, str(b1), str(0)]
                    gpu_proc = await asyncio.subprocess.create_subprocess_exec(*gpu_ecm_args,
                                                                               stdout=asyncio.subprocess.PIPE,
                                                                               stdin=asyncio.subprocess.PIPE,
                                                                               stderr=asyncio.subprocess.PIPE,
                                                                               process_group=0)
                    # synchronously write, flush, and close the one number to gpu ecm, shouldn't block
                    gpu_proc.stdin.write(str(cofactor).encode())
                    await gpu_proc.stdin.drain()
                    gpu_proc.stdin.close()


                if old_temp_save_file_path and await aios.path.isfile(old_temp_save_file_path):
                    # run cpu curves
                    residue_lines = pathlib.Path(old_temp_save_file_path).read_text().strip().split("\n")
                    cpu_ecm_args = [install_location, "-resume", "-", str(old_b1), str(old_b2)]
                    cpu_procs = []
                    input_queue = asyncio.Queue()
                    for residue_line in residue_lines:
                        await input_queue.put(residue_line)
                    tasks = []
                    for i in range(num_threads):
                        cpu_proc = await asyncio.create_subprocess_exec(*cpu_ecm_args,
                                                                        stdout=asyncio.subprocess.PIPE,
                                                                        stdin=asyncio.subprocess.PIPE,
                                                                        stderr=asyncio.subprocess.STDOUT,
                                                                        process_group=0)
                        tasks.append(asyncio.create_task(handle_stdin(input_queue, cpu_proc.stdin)))
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
                        using_lines = {}
                        for cpu_proc in cpu_procs:
                            line = await cpu_proc.stdout.readline()
                            if interrupt_level >= 3:  # stop processing CPU
                                eprint(
                                    f" {count: >{total_curves_digits}}/{total_curves}@{old_b1},{old_b2},{param}  {message if count < total_curves else '': <{len(message)+10}}")
                                cancel_tasks(tasks)
                                await kill_processes(cpu_procs)
                                await asyncio.sleep(0.2)
                                loop.stop()
                                return
                            if line:
                                line = line.decode().strip()
                                # eprint(line)
                                all_cpu_processes_done = False
                                if not line.startswith("Step 1") and not line.startswith("Input number is") and not line.startswith("Resuming ECM") and not line.startswith("GMP-ECM") and not line.startswith("Please report internal"):
                                    if line.startswith("Step 2"):
                                        count += 1
                                        eprint(f" {count: >{total_curves_digits}}/{total_curves}@{old_b1},{old_b2},{param}  {message if count < total_curves else '': <{len(message)+10}}", end="\r")
                                    elif line.startswith("Using"):
                                        using_lines[cpu_proc] = line
                                    elif line.startswith("********** Factor found in step 2:"):
                                        found_factor = int(line.strip().split(" ")[-1])
                                        if math.log10(found_factor) >= factor_size_sigma_report_threshold:
                                            eprint(f"********** BIG ECM STAGE-2 HIT: TELL YOUR FRIENDS! **********")
                                            eprint(using_lines[cpu_proc])
                                        found_factors.append(found_factor)
                                        break
                                    else:
                                        eprint(f"unknown line from CPU ECM with b1:{old_b1} b2:{old_b2}:\"{line}\"")
                                        # eprint(f"sys exit")
                                        # await shutdown(1)  # not sure how to parse yet
        await kill_gpu_procs(False)
        eprint()
        loop.stop()



if __name__ == "__main__":
    loop.create_task(main())
    loop.run_forever()
