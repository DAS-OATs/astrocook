"""Messages logged by the recipes.

Messages shared by several recipes are defined here once, so that the
same situation is always reported with the same words in the session log.

The texts taken from V1 (``astrocook.legacy.message``) were kept unchanged
with the move to V2; phase 6 reviews them with the others, by the rule of
the log messages (``ROADMAP.md``, phase 6): DEBUG for the development,
INFO for what an operation did, WARNING when the result needs a check or
something was done instead of what was asked, ERROR when there is no
result; English sentences without prefixes, values with their units. The
session log records the recipes and their parameters, not these texts.
"""

#: Logged (at the ``ERROR`` level) when a recipe cannot convert its
#: parameters, given as strings by the GUI dialogs, to the expected types;
#: the recipe then returns 0 instead of a new session.
msg_param_fail = ("I can't understand the values of the parameters. "
                  "Please try again.")
