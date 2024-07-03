import logging
import numpy as np

def parse_range(string, spec=None):
    if string=='all':
        if spec is None:
            xmin = 0
            xmax = np.infty
        else:
            xmin, xmax = np.min(spec._t['x']), np.max(spec._t['x'])
    elif string[0]=='[' and string[-1]==']':
        xmin, xmax = tuple([float(s) for s in string[1:-1].split(',')])
    else:
        logging.error("I cannot parse the range. The syntax is `[•,•]`.")
        return 0

    return xmin, xmax
