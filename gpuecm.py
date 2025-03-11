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
import time

import gmpy2
from aiofiles import os as aios

import t_level

interrupt_level = 0
gpu_proc = old_gpu_proc = cpu_procs = tasks = None
loop = asyncio.new_event_loop()
messages = ["", "", ""]


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


# tiny amount of tf to simplify ECM log searching
def tf(composite):
    cofactor = composite
    divisors = []
    for i in [2, 3, 5, 7, 11, 13]:
        while cofactor % i == 0:
            divisors.append(i)
            cofactor = cofactor // i
    return cofactor, divisors


def locate_gpu_ecm_install():
    install_location = shutil.which("ecm")
    result = subprocess.run([install_location, "-printconfig"], capture_output=True)
    compile_options = result.stdout.decode("utf-8")
    if "WITH_GPU = 1" in compile_options:
        return install_location
    eprint("Could not find GPU ECM binary location, maybe try `make install` in your ecm source dir")
    sys.exit(1)


# run a small batch to see how many curves it picks
def determine_optimal_num_gpu_curves(ecm_path, param, gpu_device):
    composite = str(797161)
    gpu_flags = ["-gpu"]
    if gpu_device:
        gpu_flags += ["-gpudevice", str(gpu_device)]
    process = subprocess.run([ecm_path] + gpu_flags + ["-param", str(param), "2", "0"], input=composite, capture_output=True, text=True)
    for line in process.stdout.split("\n"):
        match = re.search(r"(\d+) curves", line)
        if match:
            return int(match.group(1))
    assert False, "Couldn't determine optimal number of gpu curves"

def get_b1_b2_plan():
    # replace this with dynamic b1 level choice
    return [
        (145628, 2393438),
        (928607, 41925067),
        (1771522, 101780816),
        (3226570, 350917712),
        (5572142, 703952402),
        (9002399, 1450692881),
        (11381861, 2178275040),
        (17547145, 4358963387),
        (21161583, 5777034467),
        (41477251, 17337754890),
        (41477451, 17337908706),
        (49441555, 23462949664),
        (63061693, 35625232393),
        (119540944, 97480689734),
        (119040294, 96935620925),
        (118551843, 96403833440),
        (160948703, 145934182071),
        (159671931, 144434006888),
        (200832714, 192796897948),
        (334879708, 393155579274),
        (423599976, 586254865918),
        (335023594, 393468746542),
        (421896206, 582546617747),
    ]


# async task to send all the inputs, flush, and then close the stdin, i feel like this shouldn't be so difficult
async def handle_stdin(queue, stdin):
    try:
        while not queue.empty():
            to_send = (await queue.get() + "\n")
            # eprint(f"sending {to_send}", end="")
            stdin.write(to_send.encode())
            await stdin.drain()
            # queue.task_done() # not using join so I don't think we need this
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
                kill_result = proc.kill()
                if asyncio.iscoroutine(kill_result):
                    await kill_result
        except (OSError, TypeError) as e:
            if log:
                logging.exception(e)


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
    all_factors = []

    # deal with powers
    cofactor = input_number
    for i in sorted(found_factors):
        while gmpy2.is_divisible(cofactor, i):
            all_factors.append(i)
            cofactor = cofactor // i

    found_factors = all_factors

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


def t_level_lines_dict_to_lines(lines_dict):
    retlines = []
    for b1b2tup in lines_dict.keys():
        b1, b2 = b1b2tup
        curves, total_curves = lines_dict[b1b2tup]
        retlines.append([curves, b1, b2, 3])
        if curves != total_curves:
            retlines.append([total_curves - curves, b1, 0, 3])
    return retlines


async def create_gpu_proc(*args, **kwargs):
    gpu_proc = GPUProc(*args, **kwargs)
    await gpu_proc._init()
    return gpu_proc


# noinspection PyArgumentList
class GPUProc(object):

    def __init__(self, ecm_path, gpu_device, gpu_curves, save_file_path, b1, composite):
        self.ecm_path = ecm_path
        self.gpu_device = gpu_device
        self.gpu_curves = gpu_curves
        self.save_file_path = save_file_path
        self.b1 = b1
        self.composite = composite

    async def _init(self):
        gpu_device_args = ["-gpudevice", str(self.gpu_device)] if self.gpu_device else []
        gpu_curves_args = ["-gpucurves", str(self.gpu_curves)] if self.gpu_curves else []
        gpu_ecm_args = [self.ecm_path] + gpu_device_args + gpu_curves_args + ["-gpu", "-v", "-save", self.save_file_path, str(self.b1), "0"]
        self.gpu_proc = await asyncio.subprocess.create_subprocess_exec(
            *gpu_ecm_args, stdout=asyncio.subprocess.PIPE, stdin=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT, process_group=0)
        self.start_time = time.time()
        # synchronously write, flush, and close the one number to gpu ecm, shouldn't block
        self.gpu_proc.stdin.write(str(self.composite).encode())
        await self.gpu_proc.stdin.drain()
        self.gpu_proc.stdin.close()
        loop = asyncio.get_event_loop()
        self.lines = []
        self.output_parser_task = loop.create_task(self._output_parser())
        self.estimated_end_time = None

    def elapsed_time(self):
        return time.time() - self.start_time

    def remaining_time(self):
        return self.estimated_end_time - time.time() if self.estimated_end_time else None

    def estimated_total_time(self):
        return self.estimated_end_time - self.start_time  if self.estimated_end_time else None

    def has_eta(self):
        return self.estimated_end_time is not None

    def percentage_progress_string(self):
        return f"{min(100, 100 * self.elapsed_time() / self.estimated_total_time()):.01f}%" if self.estimated_end_time else ""

    def return_code(self):
        return self.gpu_proc.returncode

    async def kill(self):
        self.output_parser_task.cancel()
        return self.gpu_proc.kill()

    # Computing 11366 bits/call, 681174/10035736 (6.8%), ETA 73 + 5 = 78 seconds (~10 ms/curves)
    # Computing 11366 bits/call, 1817774/10035736 (18.1%), ETA 63 + 14 = 78 seconds (~9 ms/curves)
    # Computing 11366 bits/call, 2954374/10035736 (29.4%), ETA 55 + 23 = 77 seconds (~9 ms/curves)
    # Computing 11366 bits/call, 4090974/10035736 (40.8%), ETA 46 + 32 = 77 seconds (~9 ms/curves)
    async def _output_parser(self):
        while self.gpu_proc.returncode is None:
            line = await self.gpu_proc.stdout.readline()
            if line:
                line = line.decode()
                self.lines.append(line)
                match = re.search(r"ETA (\d+) \+ (\d+) = (\d+) seconds?", line)
                if match:
                    self.estimated_end_time = int(match.group(1)) + time.time()
            else:
                break

    async def get_output_lines(self):
        await self.output_parser_task
        return self.lines


async def main():
    global cpu_procs, gpu_proc, old_gpu_proc, tasks, interrupt_level, messages

    __version__ = "0.0.3"

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
    parser.add_argument(
        "-w",
        "--work",
        action="store",
        dest="work",
        type=float,
        help="existing t-level of work done, determines starting point of work"
    )
    parser.add_argument(
        "-p",
        "--pretest",
        action="store",
        dest="pretest",
        type=float,
        help="the desired t-level to reach in decimal digits, stops work after reaching"
    )
    parser.add_argument(
        "-t",
        "--threads",
        action="store",
        dest="threads",
        type=positive_integer,
        default=os.cpu_count(),
        help="number of threads to use in CPU stage-2"
    )
    parser.add_argument(
        "--gpudevice",
        action="store",
        dest="gpu_device",
        type=nonnegative_integer,
        help="use device <GPU_DEVICE> to execute GPU code (by default, CUDA chooses)"
    )
    parser.add_argument(
        "--gpucurves",
        action="store",
        dest="gpu_curves",
        type=positive_integer,
        help="compute on <GPU_CURVES> curves in parallel on the GPU (by default, CUDA chooses)"
    )
    parser.add_argument(
        "-o",
        "--one",
        action="store",
        dest="exit_after_one",
        help="stop ECM curves on composite after one factor is found"  # , equivalent to -x 1"
    )
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
        global interrupt_level, messages
        if signame == 'SIGINT':
             interrupt_level += 1
        if signame == 'SIGTERM':
            interrupt_level += 10
        inner_message = ""
        if interrupt_level == 1:
            inner_message = "finish GPU and CPU..."
        elif interrupt_level == 2:
            inner_message = "kill GPU, finish CPU..."
            await kill_gpu_procs()
        elif interrupt_level == 3:
            inner_message = "kill all and quit..."
        messages[2] = f"{interrupt_level}{', ' + inner_message if inner_message else ''}"
        eprint("", end="\r")  # sigint writes a newline on most terminals
        # eprint(message, end="\r")
        if interrupt_level >= 4:
            eprint("\n# Shutting Down...", end="")
            await shutdown(1)

    for signame in ('SIGINT', 'SIGTERM'):
        loop.add_signal_handler(getattr(signal, signame),
                                lambda signame=signame: asyncio.create_task(handle_signals(signame)))

    input_numbers = []
    param = 3
    factor_size_sigma_report_threshold = 60
    install_location = locate_gpu_ecm_install()
    gpu_device = args.gpu_device
    num_threads = args.threads
    exit_after_one = args.exit_after_one
    pretest = args.pretest
    work = args.work
    done_lines_dict = {}
    if work:
        done_string, _ = t_level.get_suggestion_curves_string([], 0, work, None, None, 3, 3)
        done_line = t_level.convert_string_to_parsed_lines(done_string)[0]
        done_lines_dict = {(done_line[1], done_line[2]): (done_line[0], done_line[0])}
    curves_per_batch = determine_optimal_num_gpu_curves(install_location, param, gpu_device) if not args.gpu_curves else args.gpu_curves
    with tempfile.TemporaryDirectory() as tmpdir:
        if not sys.stdin.isatty():
            input_numbers = list(map(int, sys.stdin.read().strip().split("\n")))
        gpu_proc = None
        old_gpu_proc = None
        for input_number in input_numbers:
            t_level_lines_dict = done_lines_dict
            next_lines_dict = {}
            tlev = 0
            efs = 0
            await kill_gpu_procs(True)
            gpu_proc = None
            old_gpu_proc = None
            temp_save_file_path = None
            eprint(f"# N = {input_number}", end="")
            cofactor, found_factors = tf(input_number)
            fully_factored = gmpy2.is_prime(cofactor)
            plan = get_b1_b2_plan()
            for i in range(len(plan)):
                if fully_factored or (found_factors and exit_after_one):
                    print(f"{input_number}={sorted(found_factors)}")
                    break
                if pretest and tlev >= pretest:
                    break
                b1, b2 = plan[i]
                old_b1, old_b2 = plan[i-1] if i > 0 else plan[i]
                # check if next run in plan gets us above the passed in work level, if not, skip this level
                # this is so we skip past the easy work to the productive work for our given -w <t-level>
                if work and tlev < work:
                    next_lines_dict = next_lines_dict | {(b1, b2): (curves_per_batch, curves_per_batch)}
                    next_tlev, _ = t_level.get_t_level_and_efs(t_level_lines_dict_to_lines(next_lines_dict))
                    if next_tlev < work:
                        continue
                if gpu_proc:  # wait for old gpu_proc to finish before continuing with next
                    eprint()
                    sleep_time = 0.01
                    while gpu_proc.return_code() is None:
                        tlev, efs = t_level.get_t_level_and_efs(t_level_lines_dict_to_lines(t_level_lines_dict))
                        messages[0] = f"t{tlev:0.3f}, efs: {efs:0.3f}" if tlev else ''
                        messages[1] = f"gpu: {gpu_proc.percentage_progress_string()}" if gpu_proc and gpu_proc.has_eta() else ''
                        message = f'# {", ".join(filter(bool, messages))}'
                        status_line = f"#{0: >{len(str(curves_per_batch))}}/{curves_per_batch}@{old_b1},{old_b2},{param}  {message}"
                        eprint(f"{status_line: <{shutil.get_terminal_size().columns}}", end="\r")
                        await asyncio.sleep(min(1.0, sleep_time))
                        sleep_time *= 2
                    gpu_return = gpu_proc.return_code()
                    outs = "".join(await gpu_proc.get_output_lines())
                    t_level_lines_dict[(old_b1, old_b2)] = (0, curves_per_batch)
                    tlev, efs = t_level.get_t_level_and_efs(t_level_lines_dict_to_lines(t_level_lines_dict))
                    messages[0] = f"t{tlev:0.3f}, efs: {efs:0.3f}" if tlev else ''
                    messages[1] = ''
                    message = f'# {", ".join(filter(bool, messages))}'
                    status_line = f" {0: >{len(str(curves_per_batch))}}/{curves_per_batch}@{old_b1},{old_b2},{param}  {message}"
                    eprint(f"{status_line: <{shutil.get_terminal_size().columns}}", end="\r")
                    errs = ""
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
                                eprint()
                                eprint(f"********** BIG GPU ECM STAGE-1 HIT: TELL YOUR FRIENDS! SIGMA={params[1]}:{params[2]} **********")
                                eprint(gpu_factors, end="")
                        cofactor, found_factors, fully_factored = new_factors_found(input_number, cofactor, found_factors)
                        eprint(f"# {','.join(map(str, sorted(found_factors))): <{shutil.get_terminal_size().columns}}")
                    if fully_factored or (found_factors and exit_after_one):
                        eprint()
                        print(f"{input_number}={sorted(found_factors)}")
                        break

                    # eprint(outs)
                    # eprint(errs)
                    # print(f"\n\nsigma={param}:{low_sigma}-{param}:{high_sigma}")
                    if gpu_return in [2, 6, 8]:  # found factor(s)
                        # eprint(outs)
                        # eprint(errs)
                        # don't kill here because we might need to find more factors
                        # await shutdown(1)
                        pass
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
                    gpu_proc = await create_gpu_proc(install_location, gpu_device, curves_per_batch, temp_save_file_path, b1, cofactor)

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
                    cpu_found_factors = []
                    count = 0
                    total_curves = len(residue_lines)
                    total_curves_digits = len(str(total_curves))
                    last_time = 0
                    while not all_cpu_processes_done:
                        if cpu_found_factors or (pretest and tlev >= pretest):
                            # kill cpu processes now that we have factor
                            cancel_tasks(tasks)
                            await kill_processes(cpu_procs)
                            if cpu_found_factors:
                                found_factors += cpu_found_factors
                                cofactor, found_factors, fully_factored = new_factors_found(input_number, cofactor, found_factors)
                                eprint(f"\n# {','.join(map(str, sorted(found_factors))): <{shutil.get_terminal_size().columns}}")
                            await asyncio.sleep(0.2)
                            # break out of CPU loop because factor found
                            break
                        all_cpu_processes_done = True
                        using_lines = {}
                        for cpu_proc in cpu_procs:
                            line = await cpu_proc.stdout.readline()
                            if interrupt_level >= 3:  # stop processing CPU
                                t_level_lines_dict[(old_b1, old_b2)] = (count, total_curves)
                                tlev, efs = t_level.get_t_level_and_efs(t_level_lines_dict_to_lines(t_level_lines_dict))
                                messages[0] = f"t{tlev:0.3f}, efs: {efs:0.3f}" if tlev else ''
                                messages[1] = f"gpu: {gpu_proc.percentage_progress_string()}" if gpu_proc else ''
                                message = f'# {", ".join(filter(bool, messages))}'
                                status_line = f" {count: >{len(str(curves_per_batch))}}/{curves_per_batch}@{old_b1},{old_b2},{param}  {message}"
                                eprint(f"{status_line: <{shutil.get_terminal_size().columns}}", end="\r")
                                cancel_tasks(tasks)
                                await kill_processes(cpu_procs)
                                await asyncio.sleep(0.1)
                                loop.stop()
                                return
                            if line:
                                line = line.decode().strip()
                                # eprint(line)
                                all_cpu_processes_done = False
                                if not line.startswith("Step 1") and not line.startswith("Input number is") and not line.startswith("Resuming ECM") and not line.startswith("GMP-ECM") and not line.startswith("Please report internal"):
                                    if line.startswith("Step 2"):
                                        count += 1
                                        now = time.time()
                                        if now - last_time > 0.1:
                                            t_level_lines_dict[(old_b1, old_b2)] = (count, total_curves)
                                            tlev, efs = t_level.get_t_level_and_efs(t_level_lines_dict_to_lines(t_level_lines_dict))
                                            messages[0] = f"t{tlev:0.3f}, efs: {efs:0.3f}" if tlev else ''
                                            messages[1] = f"gpu: {gpu_proc.percentage_progress_string()}" if gpu_proc and gpu_proc.has_eta() else ''
                                            message = f'# {", ".join(filter(bool, messages))}'
                                            last_time = now
                                        status_line = f" {count: >{len(str(curves_per_batch))}}/{curves_per_batch}@{old_b1},{old_b2},{param}  {message if count < total_curves else '': <{len(message) + 20}}"
                                        eprint(f"{status_line: <{shutil.get_terminal_size().columns}}", end="\r")
                                    elif line.startswith("Using"):
                                        using_lines[cpu_proc] = line
                                    elif line.startswith("********** Factor found in step 2:"):
                                        found_factor = int(line.strip().split(" ")[-1])
                                        if math.log10(found_factor) >= factor_size_sigma_report_threshold:
                                            eprint(f"********** BIG ECM STAGE-2 HIT: TELL YOUR FRIENDS! **********")
                                            eprint(using_lines[cpu_proc])
                                        cpu_found_factors.append(found_factor)
                                        found_factors.append(found_factor)
                                        break
                                    else:
                                        eprint(f"unknown line from CPU ECM with b1:{old_b1} b2:{old_b2}:\"{line}\"")
                                        # eprint(f"sys exit")
                                        # await shutdown(1)  # not sure how to parse yet
        await kill_gpu_procs(False)
        await asyncio.sleep(0.1)
        loop.stop()



if __name__ == "__main__":
    loop.create_task(main())
    loop.run_forever()
