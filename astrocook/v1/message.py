import logging
import numpy as np
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(module)s: %(message)s")
from tqdm import tqdm

msg_try_again = "Please try again."
msg_check = "Please check."

msg_exec_fail = "I can't go through with the execution."
msg_output_fail = "The output is empty. %s" % msg_try_again
msg_param_fail = "I can't understand the values of the parameters. %s" \
                  % msg_try_again
msg_param_swap = "You swapped the values of the parameters. I put them right."

def msg_attr_miss(attr):
    return "I can't find attribute %s. %s" % (attr, msg_try_again)
def msg_col_miss(col):
    return "I can't find column %s. %s" % (col, msg_try_again)
def msg_descr_miss(descr):
    return "I can't find descriptor %s. %s" % (descr, msg_try_again)
def msg_format(format):
    return "I'm importing data with %s format." % format
def msg_lim(lim):
    return "%s is badly formatted. I ignored it." % lim
def msg_empty(struct):
    return "Structure %s is empty." % struct


def msg_resol_old(check, prefix):
    if not np.logical_or(*check):
        return "[ERROR] %s: I couldn't take the resolution either from " \
               "parameter 'resol' or from the spectrum table. Please try " \
               "running Recipes > Estimate resolution before." \
               % prefix
    if check[0]:
        msg = "[INFO] %s: I've taken the resolution parameter 'resol'." \
              % prefix
        if check[1]:
            msg = msg + " Ignoring the spectrum table."
    else:
        msg = "[INFO] %s: I've taken the resolution from the spectrum table." \
              % prefix
    return msg


def msg_resol(check, prefix):
    if not np.logical_or(*check):
        return "[WARNING] %s: I couldn't take the resolution either from " \
               "parameter 'resol' or from the spectrum table. I'll estimate " \
               "it assuming 3 pixels per resolution element." \
               % prefix
    if check[0]:
        msg = "[INFO] %s: I've taken the resolution parameter 'resol'." \
              % prefix
        if check[1]:
            msg = msg + " Ignoring the spectrum table."
    else:
        msg = "[INFO] %s: I've taken the resolution from the spectrum table." \
              % prefix
    return msg



def msg_z_range(z_list):
    if len(z_list)==0:
        return "s"
    elif len(z_list)==1:
        return " at redshift %2.4f" % z_list[0]
    else:
        return "s between redshift %2.4f and %2.4f" \
                   % (np.min(z_list), np.max(z_list))

def enum_tqdm(iter, total, msg):
    return enumerate(tqdm(iter, total=total, leave=False,
                          desc="[INFO] "+msg))
