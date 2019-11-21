import logging
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(module)s: %(message)s")

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
