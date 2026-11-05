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
msg_param_fail = ("The values of the parameters cannot be read: check them "
                  "and run the recipe again.")

#: Logged (``WARNING``) when a fit runs on a spectrum without a resolution:
#: the profiles are not convolved with the instrumental profile.
msg_no_resolution = ("The spectrum has no resolution: the profiles are fitted "
                     "without the instrumental broadening, so b includes it. Set "
                     "the resolution with Edit > Set Properties.")

#: Logged (``ERROR``) when a recipe names a session that is not open; format
#: with ``name``.
msg_session_not_found = "There is no open session named '{name}'."

#: Logged (``WARNING`` or ``ERROR``) when a recipe needs the emission
#: redshift and it is not set.
msg_z_em_unset = ("The emission redshift (z_em) is not set: set it with "
                  "Edit > Set Properties.")
