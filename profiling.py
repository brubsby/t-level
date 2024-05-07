import cProfile
import importlib
import re
import gmpy2

t_level = importlib.import_module("t-level")

cProfile.run(
"""
input_string = '7557@43e6;8192@57e6,p=3;17884@11e7;8192@14e7,p=3;1532@27e7,p=3;1061@26e7;315@69e7;5415@27e7,p=3'
lines = re.split(r'(?:;|\\r?\\n)', input_string) if input_string else []
parsed_lines = list(map(t_level.parse_line, lines))
line_validations = list(map(t_level.validate_line, zip(lines, parsed_lines)))
t_level, efs = t_level.convert_lines_to_t_level_and_efs(parsed_lines)
print(t_level, efs)
""")