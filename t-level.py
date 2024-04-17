import math
import logging
import re
import argparse
import sys
import sqlite3
import os

# to build the binary, download pyinstaller with:
# pip install -U pyinstaller
# and run
# pyinstaller -F --add-data=ecmprobs.db:. t-level.py
# or
# python -m PyInstaller -F --add-data=ecmprobs.db:. t-level.py
# and find the binary in dist

conn = sqlite3.connect(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ecmprobs.db'))
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


def get_failure_probabilities(b1, curves, param):
    f = []
    if (math.isinf(b1) or math.isinf(curves)):
        return f
    c.execute(f"SELECT curves FROM ecm_probs WHERE B1 = ? AND param = ? ORDER BY curves ASC", (b1, param))
    return list(map(lambda m: pow(1.0 - (1.0/m), curves), map(lambda x: x[0], c.fetchall())))


def get_t_level(curve_b1_tuples):

    if len(curve_b1_tuples) == 0:
        return 0.0

    fp = []
    total_fp = []

    # gather failure probabilities for the given work numbers
    # make sure that the supplied b1 values are in our probability tables...
    for curves, b1, param in curve_b1_tuples:
        fp.append(get_failure_probabilities(b1, curves, param))

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
    diff = get_expected_factor_size(total_dp)
    return diff


if __name__ == "__main__":
    __license__ = "MIT"
    __version__ = "0.9.2"

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

    line_regex = r"(\d+)@(?:B1=)?(\d+e\d+|\d+)(?:,\s*(?:B2=)?(\d+e\d+|\d+))?(?:,\s*(?:(?:param|p)=)?([0-4]))?\s*"

    def parse_line(line, param=None):
        match = re.fullmatch(line_regex, line)
        if not match:
            raise ValueError(f"Malformed ecm curve string: \"{line.strip()}\"\n"
                             f"Must match {line_regex}")
        curves = sci_int(match.group(1))
        B1 = sci_int(match.group(2))
        B2 = sci_int(match.group(3))
        if param is None:
            param = sci_int(match.group(4)) if match.group(4) else 1
        logging.info(f"Curve string \"{line.strip()}\" parsed as curves={curves}, B1={B1}, B2={B2}, param={param}")
        return curves, B1, B2, param


    # TODO implement B2!=default
    def validate_line(line_tup):
        line, parsed_line = line_tup
        curves, B1, B2, param = parsed_line
        if B2:
            raise ValueError(f"Problem with curve string: \"{line.strip()}\" only ecm default B2 supported")
        return True


    def convert_lines_to_t_level(parsed_lines):
        c_at_b1_strings = list(map(lambda line: (line[0], line[1], line[3]), parsed_lines))
        return get_t_level(c_at_b1_strings)


    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"       echo <curve_string>[;<curve_string>][...] | %(prog)s [options]\n"
                    f"       printf <curve_string>[\\\\n<curve_string>][...] | %(prog)s [options]\n"
                    f"       %(prog)s [options] < <input_file>\n"
                    f"\n"
                    f"<curve_string> must full match the regex:\n"
                    f"  {line_regex}\n"
                    f"examples: 5208@11e6\n"
                    f"          5208@11e6,35133391030,1\n"
                    f"          5208@11e6,35e9,p=3\n"
                    f"          5208@B1=11e6,B2=35e9,param=0\n"
                    f"and multiple curve strings must be delimited by semicolons or newlines.")
    parser.add_argument(
        "-p",
        "--param",
        action="store",
        dest="param",
        type=int,
        help="force all input to be considered curves run using this param [0-4]"
    )
    parser.add_argument(
        "-r",
        "--precision",
        action="store",
        dest="precision",
        default=3,
        type=int,
        help="t-level decimal precision to display"
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="verbosity (-v, -vv, etc)")
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))
    args = parser.parse_args()

    if sys.stdin.isatty():
        parser.print_help()
        sys.exit(1)

    loglevel = logging.WARNING
    if args.verbose > 0:
        loglevel = logging.INFO
    if args.verbose > 1:
        loglevel = logging.DEBUG
    logging.basicConfig(level=loglevel, format="%(message)s")

    try:
        stdinput = sys.stdin.read().strip()
        lines = re.split(r'(?:;|\r?\n)', stdinput)
        parsed_lines = list(map(parse_line, lines))
        line_validations = list(map(validate_line, zip(lines, parsed_lines)))
        logging.debug(f"Validations: {line_validations}")
        t_level = convert_lines_to_t_level(parsed_lines)
        print(f"t{t_level:.{args.precision}f}")
    except ValueError as e:
        if loglevel < logging.INFO:
            logging.exception(e)
        else:
            logging.error(e)
        sys.exit(1)