import itertools
import math
import logging
import re
import argparse
import sys
import os
import sqlite3
import pathlib
import importlib.util
from . import rho
from functools import cache

# To build the binaries, run:
# python build_binaries.py


__license__ = "GPL"
__version__ = "1.0.0"

__DEFAULT_PRECISION__ = 3

if getattr(sys, 'frozen', False):
    # If the application is run as a bundle, the PyInstaller bootloader
    # extends the sys module by a flag frozen=True and sets the app 
    # path into variable _MEIPASS'.
    application_path = sys._MEIPASS
    db_path = os.path.join(application_path, 't_level', 'ecmprobs.db')
else:
    application_path = os.path.dirname(os.path.abspath(__file__))
    db_path = os.path.join(application_path, 'ecmprobs.db')

conn = sqlite3.connect(db_path, check_same_thread=False)
c = conn.cursor()


def get_differential_probability(fp_list):
    dp = []
    # calculate differential probabilities.  Numerically compute the derivative of the failure probability curve...
    dp_tmp = 0.0
    fp_len = len(fp_list)
    for i in range(len(fp_list)):
        if i == 0:
            dp.append(fp_list[1]/2.0)
        elif i == fp_len - 1:
            dp_tmp = (1.0 - fp_list[fp_len-2])/2.0
            if dp_tmp >= 0.0:
                dp.append(dp_tmp)
            else:
                dp.append(0.0)
        else:
            dp_tmp = (fp_list[i + 1] - fp_list[i - 1]) / 2.0
            if dp_tmp >= 0.0:
                dp.append(dp_tmp)
            else:
                dp.append(0.0)
    return dp


def get_expected_factor_size(dp_list):
    # calculate the expected factor size = "sum over all d-digits"( "differential probability at d-digits" * "d-digits" )
    # in our work, d-digits ranges from 10 to 100, inclusive
    diff = 0.0
    for i in range(len(dp_list)):
        # print("dp_list[" + i + "] = " + dp_list[i] + " : i+10 = " + (i+10))
        diff += dp_list[i]*(i+10)
    return diff


def get_ecm_success_probs(b1, b2, param):
    if param is None:
        param = 1
    # faster to get from cache than recalculate
    if b2 is None and (b1_level_round(b1) == b1 and 10000 <= b1 <= 50000000000):
        c.execute(f"SELECT curves FROM ecm_probs WHERE B1 = ? AND param = ? ORDER BY curves ASC", (b1, param))
        return map(lambda x: 1./x[0], c.fetchall())
    last_success_prob = None
    success_probs = []
    for digits in range(10, 101):
        if last_success_prob is not None and last_success_prob < 1 / 10000000:
            success_probs.append(0.0)
        else:
            last_success_prob = rho.ecmprob(digits, b1, b2, param)
            success_probs.append(last_success_prob)
    return success_probs


def get_failure_probabilities(b1, b2, curves, param):
    f = []
    if (math.isinf(b1) or math.isinf(curves)):
        return f
    return list(map(lambda m: pow(1.0 - m, curves), get_ecm_success_probs(b1, b2, param)))


def get_probabilities(curve_b1_tuples):
    fp = []
    total_fp = []

    # gather failure probabilities for the given work numbers
    # make sure that the supplied b1 values are in our probability tables...
    for curves, b1, b2, param in curve_b1_tuples:
        fp.append(get_failure_probabilities(b1, b2, curves, param))

    for i in range(91):
        total_fp.append(1.0)

    # combine all given failure probabilities...
    for i in range(len(fp)):
        for j in range(len(fp[i])):
            total_fp[j] = total_fp[j] * fp[i][j]

    # calculate success probabilities from the failure probabilities...
    for i in range(len(total_fp)):
        total_sp = [1.0 - fp for fp in total_fp]

    total_dp = get_differential_probability(total_fp)
    return total_fp, total_sp, total_dp


def get_t_level(curve_b1_tuples):
    return get_t_level_and_efs(curve_b1_tuples)[0]


def get_t_level_and_efs(curve_b1_tuples):

    if len(curve_b1_tuples) == 0:
        return 0.0, 0.0

    total_fp, total_sp, total_dp = get_probabilities(curve_b1_tuples)

    t_level = 0
    efs = get_expected_factor_size(total_dp)

    t_level_threshold = math.exp(-1)
    for i in range(1, len(total_fp)):
        if total_fp[i] > t_level_threshold:
            y1 = total_fp[i-1]
            y2 = total_fp[i]
            m = y2 - y1
            if m == 0:
                t_level = 0
                break
            y = t_level_threshold
            b = y2 - m * (i + 10)
            t_level = (y - b) / m
            break

    return max(0.0, float(t_level)), float(efs)


def sci_int(x):
    if x is None or type(x) in [int]:
        return x
    if type(x) != str:
        raise TypeError(f"sci_int needs a string input, instead of {type(x)} {x}")
    if x.isnumeric():
        return int(x)
    match = re.match(r"^(\d+)(?:e|[x*]10\^)(\d+)$", x)
    if not match:
        raise ValueError(f"malformed intger string {x}, could not parse into an integer")
    return int(match.group(1)) * pow(10, int(match.group(2)))


line_regex = r"\s*(\d+e\d+|\d+)(?:\/(\d+e\d+|\d+))?@(?:B1=)?(\d+e\d+|\d+)(?:,\s*(?:B2=)?(\d+e\d+|\d+))?(?:,\s*(?:(?:param|p)=)?([0-4]))?\s*"


def parse_line(line, param=None):
    if not line:
        return []
    match = re.fullmatch(line_regex, line)
    if not match:
        raise ValueError(f"Malformed ecm curve string: \"{line.strip()}\"\n"
                         f"Must match {line_regex}")
    finished_curves = sci_int(match.group(1))
    stage_1_curves = sci_int(match.group(2))
    B1 = sci_int(match.group(3))
    B2 = sci_int(match.group(4))
    if param is None:
        param = sci_int(match.group(5)) if match.group(5) else 1
    ret_curves = []
    # finished/total@B1 represents (total-finished)@B1,0 and also finished@B1
    if stage_1_curves is not None and stage_1_curves - finished_curves > 0:
        ret_curves.append((stage_1_curves - finished_curves, B1, 0, param))
    if finished_curves > 0:
        ret_curves.append((finished_curves, B1, B2, param))
    for curves, B1, B2, param in ret_curves:
        logging.info(f"Curve string: {line.strip(): <20} parsed as param={param} curves={curves: <5} B1={B1: <15} B2={B2}")
    return ret_curves


def validate_line(line_tup):
    line, parsed_line = line_tup
    curves, B1, B2, param = parsed_line
    assert(curves > 0)
    assert(B1 > 0)
    assert(B2 is None or B2 >= 0)
    assert(param in [0, 1, 2, 3, 4])
    return True


def convert_lines_to_curve_at_b1_tuples(parsed_lines):
    return list(map(lambda line: (line[0], line[1], line[2], line[3]), parsed_lines))


def convert_lines_to_t_level_and_efs(parsed_lines):
    return get_t_level_and_efs(convert_lines_to_curve_at_b1_tuples(parsed_lines))


def string_to_t_level(input_string):
    return get_t_level(convert_lines_to_curve_at_b1_tuples(convert_string_to_parsed_lines(input_string)))


def convert_string_to_parsed_lines(input_string, validation=True):
    lines = re.split(r'(?:;|\r?\n)', input_string) if input_string else []
    logging.debug(lines)
    parsed_lines = list(itertools.chain.from_iterable(map(parse_line, lines)))
    if validation:
        line_validations = list(map(validate_line, zip(lines, parsed_lines)))
        logging.debug(f"Validation results: {line_validations}")
    return parsed_lines


def convert_string_to_t_level_and_efs(input_string):
    parsed_lines = convert_string_to_parsed_lines(input_string)
    return convert_lines_to_t_level_and_efs(parsed_lines)


def get_t_level_curves(t_level, precision=__DEFAULT_PRECISION__):
    c.execute("SELECT b1, curves, MIN(ABS(curves - 10000)) FROM ecm_probs WHERE param = 1 AND digits = ?",
              (max(10, min(round(t_level), 100)),))
    b1, curves, _ = c.fetchone()
    curves = int(curves) if curves else 1
    for n in range(1, 20):
        this_t, _ = convert_lines_to_t_level_and_efs(((curves, b1, None, 1),))
        logging.debug(f"order {n} t-level estimation: {curves: >4}@{b1} = t{this_t:.{precision}f}")
        diff = t_level - this_t
        if n >= 19 or abs(diff) < pow(10, -precision)/2:
            logging.debug("precision achieved, breaking")
            return curves, b1
        last_curves = curves
        curves = max(1, int(curves * pow(2, diff / 2)))
        if curves == last_curves:
            logging.debug("stopped moving, breaking")
            return curves, b1


@cache
def get_t_level_curves_string(t_level, precision=__DEFAULT_PRECISION__):
    curves, b1 = get_t_level_curves(t_level, precision=precision)
    return f"{curves}@{b1}"


def b1_level_round(b1):
    digits = int(math.floor(math.log10(b1)))
    return int(round(b1, -digits+1))


def b1_level_string(b1):
    digits = int(math.floor(math.log10(b1)))
    return f"{b1//pow(10, digits - 1)}e{digits-1}"


@cache
def get_suggestion_curves_from_t_levels(start_t, end_t):
    return get_suggestion_curves(convert_string_to_parsed_lines(get_t_level_curves_string(start_t)), start_t, end_t, None, None, None, __DEFAULT_PRECISION__)


def get_suggestion_curves(input_lines, start_t, end_t, curves_constraint, B1_constraint, param, precision):
    # first, find some initial curves@B1 estimates for the given start_t, end_t, and constraints
    B2 = None
    if param is None:
        param = 1
    if curves_constraint is not None:
        # curves held constant, find B1 from the lookup table to get us all the way to end_t
        curves = sci_int(curves_constraint)
        c.execute("SELECT b1, MIN(ABS(curves - ?)) FROM ecm_probs WHERE param = 1 AND digits = ?",
                  (curves, max(10, min(round(end_t), 100))))
        B1, _ = c.fetchone()
    else:
        if B1_constraint is not None:
            # B1 held constant (strange but ok)
            B1 = b1_level_round(sci_int(B1_constraint))
        else:
            # no constraints, choose B1 first from the regression formula
            B1 = max(b1_level_round(get_regression_b1_for_t(end_t)), 10000)
        diff_t = end_t - start_t
        # diff = 2 is the point at which the work required will be roughly double
        if diff_t >= 2:
            # if we're doing more than twice the amount of work already done, the lookup table will
            # be closer to the end curve amount, as it holds the curves to get from t0 to end_t

            # make sure we get a B1 that is in our lookup tables
            c.execute(f"SELECT MAX(B1) FROM ecm_probs WHERE B1 <= ? AND param = ?", (B1, param))
            B1 = c.fetchone()[0]
            lookup_t = max(10, min(math.floor(end_t), 100))
            c.execute("SELECT curves FROM ecm_probs WHERE param = ? AND digits = ? AND B1 = ?",
                      (param, lookup_t, B1))
            curves = round(c.fetchone()[0])

            # if start_t is in the lookup table, subtract those curves from our initial estimate curves
            if 10 <= start_t <= 100:
                c.execute("SELECT curves FROM ecm_probs WHERE param = ? AND digits = ? AND B1 = ?",
                          (param, round(start_t), B1))
                curves = max(1, round(curves - c.fetchone()[0]))
        else:
            # for diff_t < 2, the incremental table at https://members.loria.fr/PZimmermann/records/ecm/params.html
            # will provide roughly the amount of curves to do log10(2^5) = ~1.505 t-levels of work.
            curves = max(1, round(diff_t * pow(B1 / 150, 2 / 3)))

    # next, improve upon the initial guess iteratively
    last_diff = last_t_level = last_curves = last_B1 = diff = t_level = None
    start_curves = curves
    for n in range(0, 20):
        trial_lines = [*input_lines, (curves, B1, B2, param)]
        if t_level != None:
            last_t_level = t_level
        if diff != None:
            last_diff = diff
        t_level, _ = convert_lines_to_t_level_and_efs(trial_lines)
        diff = end_t - t_level
        if last_diff and abs(last_diff) < abs(diff):
            curves = last_curves
            B1 = last_B1
            diff = last_diff
            t_level = last_t_level
        else:
            logging.debug(
                f"order {n: >2} t-level estimation: {curves: >4}@{b1_level_string(B1)} = t{t_level:.{precision}f}, diff={diff}")
        if abs(diff) < pow(10, -precision - 1):
            logging.debug("precision achieved, breaking")
            break
        last_curves = curves
        last_B1 = B1
        if curves_constraint:
            # adjust B1 level
            mult = 1 + pow(2, -n+5)
            B1 *= pow(mult, (1 if diff > 0 else -1))
            B1 = min(max(1000, b1_level_round(B1)), int(50e9))
            if B1 == last_B1 and mult != 1:
                break
        else:
            # adjust curves first, maybe then B1 if no constraint
            addend = start_curves * pow(2, -n-1)
            curves += addend * (1 if diff > 0 else -1)
            curves = max(1, round(curves))
            if curves == last_curves and addend != 0:
                logging.debug(f"terminal curves achieved on iteration {n}: {curves}")
                break
    return curves, B1, B2, param, t_level


def get_suggestion_curves_string(input_lines, start_t, end_t, curves_constraint, B1_constraint, param, precision):
    curves, B1, B2, param, t_level = get_suggestion_curves(input_lines, start_t, end_t, curves_constraint, B1_constraint, param, precision)
    B2 = "" if B2 is None else f",{B2}"
    p = "" if param == 1 else f",p={param}"
    return f"{curves}@{b1_level_string(B1)}{B2}{p}", t_level


def get_multiple_suggestion_curves_string(input_lines, start_t, end_t, interval_t, curves_constraint, param, precision):
    ret_strs = []
    t_level = start_t
    num_steps = math.ceil(round((end_t - start_t) / interval_t, 1))
    t_range = [start_t + i * interval_t for i in range(num_steps)]
    for inner_start_t in t_range:
        inner_end_t = min(end_t, inner_start_t + interval_t)
        curves, B1, B2, param, t_level = get_suggestion_curves(input_lines, inner_start_t, inner_end_t, curves_constraint, None, param, precision)
        input_lines.append((curves, B1, B2, param))
        B2 = "" if B2 is None else f",{B2}"
        p = "" if param == 1 else f",p={param}"
        ret_strs.append(f"{curves}@{b1_level_string(B1)}{B2}{p}")
    return ";".join(ret_strs), t_level


def get_regression_b1_for_t(end_t):
    # from https://members.loria.fr/PZimmermann/records/ecm/params.html
    B1 = int(math.exp(0.0750 * math.log2(pow(10, end_t)) + 5.332))
    return B1


def graph_parsed_lines(work_lines, parsed_lines):
    work_fp, work_sp, work_dp = get_probabilities(convert_lines_to_curve_at_b1_tuples(work_lines))
    next_fp, next_sp, next_dp = get_probabilities(convert_lines_to_curve_at_b1_tuples(parsed_lines))
    x = list(range(10,101))
    digits = 100
    composite = pow(10,digits)
    rho.rhoinit(256,30)
    dickman_mu = [rho.dickmanmu(digits/i, digits/i , composite) for i in x]
    mu_prior = [dickman_mu[i] * work_fp[i] for i in range(len(next_sp))]
    merten = [0.56146 for i in x]
    dickman_local = [rho.dickmanrho(digits/i) for i in x]
    work_dickman_local = [dickman_local[i] * work_fp[i] for i in range(len(next_sp))]
    next_dickman_local = [dickman_local[i] * next_fp[i] for i in range(len(next_sp))]
    diff_dickman_local = [work_dickman_local[i] - next_dickman_local[i] for i in range(len(next_sp))]
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams["figure.figsize"] = (10, 3)
    plt.xticks(range(10, 101, 5))
    plt.yticks(list(map(lambda x: x * 0.1, range(0, 11, 1))))
    plt.plot(x, work_fp, label="prev fail prob")
    plt.plot(x, work_sp, label="prev success prob")
    # plt.plot(x, work_dp, label="work dp")
    plt.plot(x, next_fp, label="next failure prob")
    plt.plot(x, next_sp, label="next success prob")
    # plt.plot(x, next_dp, label="dp")
    plt.plot(x, dickman_mu, label="rho mu")
    plt.plot(x, dickman_local, label=f"rho")
    # plt.plot(x, work_dickman_local, label="workdick")
    # plt.plot(x, next_dickman_local, label="nextdick")
    # plt.plot(x, diff_dickman_local, label="diffdick")
    # plt.plot(x, mu_prior, label="mu prior")
    plt.legend()
    plt.tight_layout()
    plt.xlim(10, 100)
    plt.ylim(0, 1)
    plt.xlabel("p-digits")
    plt.ylabel("probability")
    # plt.savefig('B2_timing.png')
    plt.show()
    pass


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"       echo <curve_string>[;<curve_string>][...] | %(prog)s [options]\n"
                    f"       printf <curve_string>[\\\\n<curve_string>][...] | %(prog)s [options]\n"
                    f"       %(prog)s [options] < <input_file>\n"
                    f"       %(prog)s [options] -i <input_file>\n"
                    f"       %(prog)s [options] -q\"<curve_string>[;<curve_string>][...]\"\n"
                    f"\n"
                    f"<curve_string> must full match the regex:\n"
                    f"  {line_regex}\n"
                    f"examples: 5208@11e6\n"
                    f"          5208@11e6,35133391030,1\n"
                    f"          5208@11e6,35e9,p=3\n"
                    f"          5208@B1=11e6,B2=35e9,param=0\n"
                    f"          5208/8192@11e6  # 5208 finished curves, 2984 stage-1 curves \n"
                    f"and multiple curve strings must be delimited by semicolons or newlines.")
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
    parser.add_argument(
        "-q",
        type=str,
        action="store",
        dest="expression",
        help="direct curve strings expression input",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        action="store",
        dest="filename",
        help="curve strings file input",
    )
    parser.add_argument(
        "-w",
        "--work",
        action="store",
        dest="work",
        type=float,
        help="existing t-level of work done, to add input curves to"
    )
    # parser.add_argument(
    #     "-f",
    #     "--force-param",
    #     action="store",
    #     dest="force_param",
    #     type=int,
    #     help="force all input to be considered curves run using this ecm param [0-4]"
    # )
    parser.add_argument(
        "-r",
        "--precision",
        action="store",
        dest="precision",
        default=__DEFAULT_PRECISION__,
        type=int,
        help="t-level decimal precision to display"
    )
    parser.add_argument(
        "-e",
        "--efs",
        action="store_true",
        dest="efs",
        help="also display expected factor size calculation",
    )
    parser.add_argument(
        "-t",
        action="store",
        dest="t_level",
        type=float,
        help="the desired t-level to reach, suggests curve string to run to achieve the given t-level"
    )
    parser.add_argument(
        "-n",
        "--work-interval",
        action="store",
        dest="work_interval",
        type=float,
        help="when suggesting curve strings to run to achieve the provided -t level, split them up into sections of "
             "t-level work this big"
    )
    parser.add_argument(
        "-c",
        "--curves",
        action="store",
        dest="curves",
        type=int,
        help="constrain suggestion curve quantity to a value, use with -t"
    )
    parser.add_argument(
        "-p",
        "--param",
        action="store",
        dest="param",
        type=int,
        help="constrain suggested param value to a value [0-4], default 1, use with -t"
    )
    parser.add_argument(
        "-b",
        action="store",
        dest="B1",
        type=str,
        help="constrain suggested B1 to a value, use with -t"
    )

    graphing_available = (importlib.util.find_spec("matplotlib") is not None and
                          importlib.util.find_spec("numpy") is not None)

    if graphing_available:
        parser.add_argument(
            "-g",
            "--graph",
            action="store_true",
            dest="graph",
            help="use matplotlib to visualize the probability density functions"
        )

    args = parser.parse_args()

    loglevel = logging.WARNING
    if args.verbose > 0:
        loglevel = logging.INFO
    if args.verbose > 1:
        loglevel = logging.DEBUG
    logging.basicConfig(level=loglevel, format="%(message)s")

    curve_inputs = []
    if args.expression:
        curve_inputs.append(args.expression)
    if not sys.stdin.isatty():
        curve_inputs.append(sys.stdin.read().strip())
    if args.filename:
        try:
            file_input = pathlib.Path(args.filename).read_text().strip()
            # logging.debug(file_input)
            curve_inputs.append(file_input)
        except FileNotFoundError as e:
            logging.error(e)
            sys.exit(1)
    if not curve_inputs and not args.t_level and sys.stdin.isatty():
        parser.print_help()
        sys.exit(1)

    if args.t_level is None:
        if args.work_interval is not None or args.B1 is not None or args.param is not None or args.curves is not None:
            print("-n, -c, -b, and -p flags can only be used with -t flag, exiting...")
            sys.exit(1)
        if args.B1 is not None and args.curve is not None:
            print("-c and -b flags cannot both be used to constrain both curves and B1 simultaneously, exiting...")
            sys.exit(1)
        if args.B1 is not None and args.work_interval is not None :
            print("-n and -b flags cannot both be used to constrain both work interval and B1 simultaneously, exiting...")
            sys.exit(1)

    if args.work:
        if args.work < 5 or args.work > 100:
            print("-w flag must be >= 5 and <= 100")
            sys.exit(1)
        work_input = get_t_level_curves_string(args.work, args.precision)
        curve_inputs.append(work_input)

    input_string = "\n".join(curve_inputs).strip()
    try:
        parsed_lines = convert_string_to_parsed_lines(input_string)
        if getattr(args, 'graph', False):
            parsed_work_lines = convert_string_to_parsed_lines(work_input) if args.work else ""
            graph_parsed_lines(parsed_work_lines, parsed_lines)
        t_level, efs = convert_lines_to_t_level_and_efs(parsed_lines)
    except ValueError as e:
        if loglevel < logging.INFO:
            logging.exception(e)
        else:
            logging.error(e)
        sys.exit(1)
    if args.expression or not sys.stdin.isatty() or args.filename:
        print(f"t{t_level:.{args.precision}f}")
        if args.efs:
            print(f"efs:{efs:.{args.precision}f}")
    if args.t_level:
        if args.work_interval:
            curves, new_t_level = \
                get_multiple_suggestion_curves_string(parsed_lines, t_level, args.t_level, args.work_interval,
                                                      args.curves, args.param, args.precision)
        else:
            curves, new_t_level = \
                get_suggestion_curves_string(parsed_lines, t_level, args.t_level, args.curves,
                                             args.B1, args.param, args.precision)
        print(f"Running the following will get you to t{new_t_level:.{args.precision}f}:")
        print(f"{curves}")


if __name__ == "__main__":
    main()
