import logging
import numpy as np
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(module)s: %(message)s")
from tqdm import tqdm

msg_exec_fail = "I can't go through with the execution."
msg_output_fail = "The output is empty. Please try again."
msg_param_fail = "I can't understand the values of the parameters. "\
                 "Please try again."
msg_param_swap = "You swapped the values of the parameters. I put them right."

def msg_attr_miss(attr):
    return 'Attribute %s is missing.' % attr
def msg_descr_miss(descr):
    return 'Descriptor %s is missing.' % descr
def msg_format(format):
    return "I'm importing data with %s format." % format

def msg_z_range(z_list):
    if len(z_list)==0:
        return "s"
    elif len(z_list)==1:
        return " at redshift %2.4f" % z_list[0]
    else:
        return "s between redshift %2.4f and %2.4f" \
                   % (np.min(z_list), np.max(z_list))

def enum_tqdm(iter, total, msg):
    return enumerate(tqdm(iter, ncols=120, total=total, leave=False,
                          desc="[INFO] "+msg))
