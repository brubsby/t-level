import math
import re
import subprocess
import t_level
import re

import numpy as np

from scipy import interpolate
from scipy.optimize import bisect, curve_fit, minimize

B1_timing = None
B1_GPU_timing = None
B2_timing = None


def load_B1_timing(composite):
    """Load CPU Step 1 timing"""
    # take as command line option
    global B1_timing
    if B1_timing is None:
        B1 = 100000
        if type(composite) != bytes:
            composite = bytes(str(composite), "utf-8")
        process = subprocess.run(["ecm", "-param", "3", str(B1), "0"], input=composite,
                                 capture_output=True)
        output = process.stdout.split(b"\n")
        timing = float(re.match(".*?([0-9]+)ms", str(output[-2], 'utf-8')).group(1))
        B1_timing = timing/B1/1000
        print(f"B1 CPU Timing: {timing}")
    return B1_timing


def load_GPU_B1_timing(composite, gpucurves):
    """Load GPU Step 1 timing"""
    # take as command line option
    global B1_GPU_timing
    if B1_GPU_timing is None:
        B1 = 100000
        if type(composite) != bytes:
            composite = bytes(str(composite), "utf-8")
        process = subprocess.run(["ecm", "-gpu", "-gpucurves", str(gpucurves), str(B1), "0"], input=composite,
                                 capture_output=True)
        output = process.stdout.split(b"\n")
        timing = float(re.match(".*?([0-9]+)ms of GPU time", str(output[-2], 'utf-8')).group(1))
        B1_GPU_timing = timing/B1/1000/gpucurves
        print(f"B1 GPU Timing: {timing}")
    return B1_GPU_timing


def load_B2_timing():
    """Load CPU Step 2 timings at various B2 values"""

    # (B2, time(ms))
    B2_timing = []

    # TODO add more values to the log or something
    with open("B2_timing.log") as f:
        B2 = None
        for line in f:
            match = re.match("Using B1=[0-9]*, B2=([1-9][0-9]*), ", line)
            if match:
                B2 = int(match.group(1))
                #print (B2, "\t", line.strip())

            match = re.match("Step 2 took ([0-9]+)ms", line)
            if match:
                assert B2
                timing = int(match.group(1)) / 1000
                #print (timing, B2, "\t", line.strip())
                B2_timing.append((B2, timing))
                B2 = None

    return sorted(B2_timing)


def B2_timing_guess(timings):
    """Take an educated guess at B2 timing based on interpolation between values"""

    x, y = tuple(zip(*timings))
    assert x == tuple(sorted(x))
    # Handle in-domain estimates
    f = interpolate.interp1d(x, y)

    # Handles out-of-domain estimates
    def powerlaw(x, c, e): return c * x ** e

    OUT_START = len(x)-4
    power_params = curve_fit(powerlaw, x[OUT_START:], y[OUT_START:])[0]
    g = lambda x: powerlaw(x, *power_params)

    def h(b2):
        return f(b2) if b2 < max(x) else g(b2)

    '''
    print(f"X from {x[1]:.1e} to {x[-1]:.1e}")
    print(f"Y from {y[1]} to {y[-1]}")

    import matplotlib.pyplot as plt
    import numpy as np
    new_x = np.logspace(np.log10(x[0]), np.log10(x[-1]), num=200, endpoint=False)
    ext_x = np.logspace(np.log10(x[OUT_START]), np.log10(max(x))+2, num=200)

    plt.plot(x, y, 'x', label="measured")
    plt.plot(new_x, f(new_x), label="interpolation")
    plt.plot(ext_x, g(ext_x), label="estimate (out)")
    plt.xlim(new_x[0], ext_x[-1])
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.tight_layout()
    plt.savefig('B2_timing.png')
    plt.show()
    # '''

    return h


def estimate_remaining_curve_proportion(done_t_level, target_t_level):
    assert done_t_level <= target_t_level
    d = pow(2, (target_t_level - done_t_level) / 2)
    return (d - 1) / d


def add_t_levels(t1, t2, precision=3):
    t, _ = t_level.convert_string_to_t_level_and_efs(
        t_level.get_t_level_curves(t1, precision) + ";" + t_level.get_t_level_curves(t2, precision))
    return t


def get_remaining_curves_to_t_level(done_t_level, target_t_level, total_curves):
    return math.ceil(estimate_remaining_curve_proportion(done_t_level, target_t_level) * total_curves)


def estimate_t_level_after_run(c, B1, B2, p, prev_curves_string):
    if p == None:
        p = 3
    curvestring = f"{c}@{B1},{B2},{p};{prev_curves_string}"
    t, efs = t_level.convert_string_to_t_level_and_efs(curvestring)
    return t


def estimate_t_diff_for_run(c, B1, B2, p, prev_curves_string, existing_work=None):
    if existing_work is None:
        existing_work, _ = t_level.convert_string_to_t_level_and_efs(prev_curves_string)
    return estimate_t_level_after_run(c, B1, B2, p, prev_curves_string) - existing_work


def optimize_t(previous_curve_strings, GPU_SPEEDUP, CPU_CORES, PARAMETER, composite):
    B1_timing = load_B1_timing(composite)
    B2_timings = load_B2_timing()

    B1_time_func = lambda B1: B1 * B1_timing
    B2_time_func = B2_timing_guess(B2_timings)

    def time_for_tX(B1, B2):
        curves = 8192
        B1_time = B1_time_func(B1) / GPU_SPEEDUP
        B2_time = (B2_time_func(B2) - B2_time_func(B1)) / CPU_CORES
        return curves, curves * (B1_time + B2_time), curves * max(B1_time, B2_time)

    def constrain_B2(B1):
        """Find B2 that takes B1 / B1_speedup time"""
        B1_time = B1_time_func(B1) / GPU_SPEEDUP
        B2_time_goal = B1_time
        #print(f"{B1:<10}\t{B1_time:0.5f}\t{B2_time_goal:0.5f}\t{B2_time_func(B1):0.5f}\t{B2_time_func(20000*B1):0.5f}")
        def test(B2):
            return (B2_time_func(B2) - B2_time_func(B1)) / CPU_CORES - B2_time_goal

        # B2 = int(minimize(test, B1, bounds=[[B1, np.inf]], method='SLSQP').x.item())
        # return B2
        MAX_B1 = 10 ** 20
        if test(MAX_B1) < 0:
            return MAX_B1

        t = bisect(test, B1, MAX_B1, maxiter=1000)
        return int(t)


    # initial guess based off of normal B1 rules
    starting_t_level, efs = t_level.convert_string_to_t_level_and_efs(previous_curve_strings)
    B1_best = max(1000000, int(179.2694750411 * pow(1.2684228593, efs)))
    B2_best = constrain_B2(B1_best)
    curves_best, _, time_best = time_for_tX(B1_best, B2_best)
    t_diff_best = estimate_t_diff_for_run(curves_best, B1_best, B2_best, PARAMETER, previous_curve_strings, existing_work=starting_t_level)
    t_diff_per_time_best = t_diff_best/time_best

    # This is way less efficient than binary search but easier to code
    for i in range(100):
        one_pct = B1_best // 100
        for B1_test in (one_pct * 150, one_pct * 125, one_pct * 113, one_pct * 103, int(one_pct*100.5), int(one_pct*100.05),
               one_pct * 98, one_pct * 99, int(one_pct*99.9), int(one_pct*99.99)):
            B2_test = constrain_B2(B1_test)

            curves_test, _, time_test = time_for_tX(B1_test, B2_test)
            t_diff_test = estimate_t_diff_for_run(curves_test, B1_test, B2_test, PARAMETER, previous_curve_strings, existing_work=starting_t_level)
            t_diff_per_time_test = t_diff_test/time_test

            if t_diff_per_time_test > t_diff_per_time_best:
                B1_best = B1_test
                B2_best = B2_test
                time_best = time_test
                t_diff_best = t_diff_test
                t_diff_per_time_best = t_diff_per_time_test
                # print(f"Curves={curves_best} B1={B1_best:<15} B2={B2_best:<15} Time={time_best:0.5f}s t_diff={t_diff_test:0.5f} t/s={t_diff_per_time_test:0.9f}")
                break
        else:
            break


    B2_ratio = B2_best / B1_best
    curves, _, _ = time_for_tX(B1_best, B2_best)
    B1_time = curves * B1_time_func(B1_best) / GPU_SPEEDUP
    B2_time = curves * (B2_time_func(B2_best) - B2_time_func(B1_best)) / CPU_CORES
    stg_1_t_level, _ = t_level.convert_string_to_t_level_and_efs(f"0/{curves}@{B1_best},{B2_best},{PARAMETER}{';'+previous_curve_strings if previous_curve_strings else ''}")
    previous_curve_strings = f"{previous_curve_strings+';' if previous_curve_strings else ''}{curves}@{B1_best},{B2_best},{PARAMETER}"
    total_t_level, _ = t_level.convert_string_to_t_level_and_efs(previous_curve_strings)
    if False:
        print("B1={}, B2={}   {} curves".format(B1_best, B2_best, curves))
        print("Step 1 takes: {:0.5f} seconds".format(B1_time))
        print("Step 2 takes: {:0.5f} seconds".format(B2_time))
        print("Running on GPU + {} cores: {:0.5f} seconds".format(CPU_CORES, max(B1_time, B2_time)))
        print("Delta t-level is {:0.5f}".format(t_diff_best))
        print("Total t-level is {:0.5f}".format(total_t_level))
    else:
        print ("| {:<0.4f} | {:<0.4f} | {:13} | {:19} | {:4.0f} | {:6} | {:12.1f}s | {:12.1f}s |".format(
            stg_1_t_level, total_t_level, B1_best, B2_best, B2_ratio, curves, B1_time, B2_time), end="")

    return B1_best, B2_best, total_t_level, max(B1_time, B2_time), previous_curve_strings


def optimize():
    COMPOSITE = 135512287471878204193633596210376549924578479626664178122249113294162659711944754598503615020622210846798939754751372353713360604667619430324232979216887721550925726303
    GPU_CORES = 8192
    CPU_CORES = 16
    PARAMETER = 3
    STARTINGT = 50
    MAXIMUM_T = 65
    GPU_SPEEDUP = load_B1_timing(COMPOSITE) / load_GPU_B1_timing(COMPOSITE, GPU_CORES)
    print(f"GPU Speedup: {GPU_SPEEDUP:.0f}")
    UPPER_SEP = "|-------------------|------------------------------------------------------------------------------------------------------|"
    LOWER_SEP = "|---------|---------|---------------|---------------------|------|--------|---------------|---------------|----------------|"
    print(UPPER_SEP)
    print(      "|  cumulative work  |                                                                                                      |")
    print(LOWER_SEP)
    print(      "| stg-1-t | stg-2-t |    optimal B1 |          optimal B2 |Bratio| curves |  stage-1 time |  stage-2 time |     total time |")
    print(LOWER_SEP)
    t_total = STARTINGT if STARTINGT > 15 else 0
    sum_time = 0
    prev_curve_strings = "" if t_total == 0 else t_level.get_t_level_curves(t_total, 3)
    if t_total > 0:
        curves, B1 = prev_curve_strings.split("@")
        print("| {:<0.4f} | {:<0.4f} | {:13} | {:19} | {:4} | {:6} | {:12.1f}s | {:12.1f}s | {:13.1f}s |".format(
            t_total, t_total, int(B1), "", "", int(curves), 0, 0, 0))
    while t_total < MAXIMUM_T:
        _, _, t_total, rowtime, prev_curve_strings = optimize_t(prev_curve_strings, GPU_SPEEDUP, CPU_CORES, PARAMETER, COMPOSITE)
        sum_time += rowtime
        print(f" {sum_time:13.1f}s |")

    print(LOWER_SEP)
    print()
    print(prev_curve_strings)
    print()
    for curve_string in prev_curve_strings.split(';'):
        split = re.split(r'[,@]', curve_string)
        if len(split) > 2:
            print(split[1], split[2])

optimize()
# ests = []
# for i in [2303    , 21558   , 235581  , 2947256 ]:
#     ests.append(estimate_t_level_for_run(i, 1533946, 18097830,  existing_work=0))
# print(ests)