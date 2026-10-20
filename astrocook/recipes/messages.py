"""Messages logged by the recipes.

Messages shared by several recipes are defined here once, so that the
same situation is always reported with the same words in the session log.

The texts taken from V1 (``astrocook.legacy.message``) are unchanged, so
that the log of a session does not change with the move to V2.
"""

#: Logged (at the ``ERROR`` level) when a recipe cannot convert its
#: parameters, given as strings by the GUI dialogs, to the expected types;
#: the recipe then returns 0 instead of a new session.
msg_param_fail = ("I can't understand the values of the parameters. "
                  "Please try again.")
