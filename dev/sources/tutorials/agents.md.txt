# Using Astrocook from an AI Agent

:::{warning}
The server is **being completed**: since the pilot every recipe is a
tool, and tools to inspect and edit the sessions are coming. Tools,
parameters and outputs may change in later versions.
:::

An AI agent, such as [Claude Code](https://claude.com/claude-code) or the
Claude desktop app, can drive Astrocook for you: you ask in plain words
("estimate the continuum and fit the CIV doublets"), the agent calls the
Astrocook recipes, reads their results and tells you what it found. You
then review its work in the GUI.

The agent talks to Astrocook through the **Model Context Protocol (MCP)**,
an open standard for connecting AI applications to other programs.
Astrocook provides an MCP *server*: a small program, `astrocook-mcp`,
that the agent (the *client*) starts on your computer and that offers it
a list of *tools*, each a recipe with a description and typed
parameters. Everything runs locally: your data are read from your disk
and stay there, though what the tools return (summaries, not the
spectra) is sent to the AI model, as anything else in the conversation.

## 1. Installation

The server is an optional part of Astrocook, installed with the `[mcp]`
extra (Python 3.12 or later, preferably in a virtual environment, see
[Installation](../installation)):
```
pip install "astrocook[mcp] @ git+https://github.com/DAS-OATs/astrocook.git"
```
This also installs the `mcp` package (the official Python SDK of the
protocol) and the `astrocook-mcp` command. The macOS app does not include
the server.

Find the full path of the command, which the client needs (it does not
see your virtual environment):
```
which astrocook-mcp     # On Windows: where astrocook-mcp
```
e.g. `/Users/me/ac2/bin/astrocook-mcp`.

## 2. Registering the Server

**Claude Code** (the `claude` command; if it is not installed, see
[its setup page](https://code.claude.com/docs/en/setup); on macOS and
Linux: `curl -fsSL https://claude.ai/install.sh | bash`):
```
claude mcp add astrocook -- /Users/me/ac2/bin/astrocook-mcp
```
To have it in all your projects, add `--scope user` before the name
(`claude mcp add --scope user astrocook -- ...`). `claude mcp list`
shows the registered servers; in a session, `/mcp` shows their state.

**Claude desktop app:** add the server to its configuration file
(Settings > Developer > Edit Config; on macOS
`~/Library/Application Support/Claude/claude_desktop_config.json`, on
Windows `%APPDATA%\Claude\claude_desktop_config.json`), then restart the
app:
```json
{
  "mcpServers": {
    "astrocook": {
      "command": "/Users/me/ac2/bin/astrocook-mcp"
    }
  }
}
```

**Other clients:** any MCP client that can start a local server over
*stdio* works the same way: the command is `astrocook-mcp`, with no
arguments.

## 3. The Tools

**Every recipe is a tool**, named after the recipe, prefixed with its
menu: *Continuum > Auto-estimate Continuum* is `continuum_estimate_auto`,
*Absorbers > Auto-detect Doublets* is `absorbers_doublets_auto`, and so on
for the menus *File*, *Edit*, *Flux*, *Continuum* and *Absorbers*, also the
recipes that the GUI runs from the System Inspector (e.g.
`absorbers_fit_component`). The tools are made from the same description
of the recipes as the dialogs of the GUI: the same parameters, defaults,
bounds and texts, and a recipe added to Astrocook becomes a tool by
itself. Other tools handle the sessions:

| Tool | What it does |
|---|---|
| `session_open` | Opens a spectrum or a session file (`.acs2`, V1 `.acs`) as a new session; optionally sets z_em and the resolving power (*Edit > Set Properties*). |
| `session_save` | Saves the session as `.acs2`, with the original data and the session log. |
| `session_list` | Lists the open sessions: where each comes from, its spectrum in a few numbers, its history, whether it has changes not saved. |
| `session_describe` | Describes a session: the columns of the spectrum, the components per series, the log that came with the file. |
| `session_history` | Lists the states of a session: as opened, then one per recipe, with its parameters and what it changed. |
| `session_undo`, `session_redo` | Go back to the state before the last recipe (or several), and forward again, as *Edit > Undo* and *Redo* in the GUI. |
| `session_duplicate` | Makes a copy of a session, with its history and log, to try something without touching the original. |
| `session_close` | Closes a session; refuses if it has changes not saved, unless asked to discard them. |

**Nothing is lost.** Every state of a session is kept, as in the GUI: if
a recipe removes components (e.g. *Auto-detect Doublets* replaces those
on the red side), its result says so, with the reason for each, and
`session_undo` brings them back without running anything again. The log
follows the state: a session saved after an undo replays without the
recipe undone.

Every recipe tool takes the `session_id` of the session to run on, the
parameters of the recipe, and:
- `provide_missing`: see "What a recipe needs" below;
- `max_listed`: the maximum number of items in each list of the output.

It returns the session after the call (the spectrum in a few numbers, the
number of log entries), what the recipe did in one line (`summary`) and
what else is worth knowing (`notes`), the warnings, and what changed,
listed: the components removed (with the reason, when the recipe gives
one) and changed, the columns and properties of the session. Some results
add a summary of their own: the continuum and the absorbed pixels after a
recipe that makes a continuum, the regions and candidate systems after an
identification, the new components (with errors, reduced chi2 and flags)
after a fit.

A recipe that opens new sessions in the GUI (e.g. *Edit > Split*, or
*Absorbers > Identify Absorption Lines*) does so here too: the session it
runs on does not change, and the result lists the new ones, each with its
`session_id`, named after the original (`J0103-1305-1`, ...). A parameter
that names another open session (e.g. the reference of *Flux > Equalize*)
takes its `session_id`.

*File > Export to ASCII* (`file_export_ascii`) writes the files in the
folder `path`, named after `prefix`, and does not replace files that exist
unless asked (`overwrite=true`).

The sessions live in the server while the client is running: save the
result with `session_save` before closing the client.

## 4. An Example

A request like this one:

> Open `/data/J0103-1305.acs2`, with z_em = 4.072 and R = 45000. Estimate
> the continuum, fit the metal doublets on the red side, and save the
> result as `/data/J0103-1305_agent.acs2`. Which components are not
> reliable?

leads the agent to call `session_open`, `continuum_estimate_auto`,
`absorbers_doublets_auto` and `session_save` in turn. Each call returns
a structured summary, for instance (from `continuum_estimate_auto`):
```json
"continuum": {
  "absorbed_fraction": 0.648,
  "absorbed_fraction_forest": 0.862,
  "absorbed_fraction_red": 0.207,
  "snr_median": 29.6
}
```
and, from `absorbers_doublets_auto`, the number of components added and
removed and the list of the new ones, each with its flags (e.g.
`b_at_lower_bound`), from which the agent answers the last question.

Each result also says what the recipe did in one line (`summary`) and
what else is worth knowing (`notes`): the components removed and why,
those whose fit failed, the steps run first.

If something goes wrong (a missing file, a value that cannot be used), the
tool returns Astrocook's error message, which names the recipe and the
reason, and the agent can usually correct its request by itself.

**What a recipe needs.** A recipe that needs something the session does
not have (the emission redshift, the resolution, a continuum) does not
run: the tool says what is missing and which recipe provides it, and the
agent decides. z_em and the resolution are given to `session_open` or
`edit_set_properties`. A
continuum is made with `continuum_estimate_auto`, or, if the agent calls
the tool again with `provide_missing=true`, with the default parameters of
*Auto-estimate Continuum* before the recipe: the step is listed in
`provided` and recorded in the session log, before the recipe, as if it
had been run on its own.

## 5. Reviewing the Result in the GUI

Open the saved archive in the GUI (`astrocook /data/J0103-1305_agent.acs2`).
The session log (*View > View Session Log*) lists the recipes the agent ran, with
their parameters, as if you had run them from the menus: you can inspect
the components in the System Inspector, correct them, refit them, or
replay the log as a script.

## 6. Limits of the Pilot

- No tools yet to list the components or read the spectrum (coming).
- Paths are on the computer where the server runs; give them in full.
- `absorbers_doublets_auto` replaces the components already on the red
  side. It finds the absorbed regions on a new mask made with the default
  parameters (`mask="new"`, used only inside), or on the mask of the
  session (`mask="session"`, e.g. that of `continuum_estimate_auto` with
  other parameters). A call can take minutes on a large spectrum.
- The recipes are the same as in the GUI, with the same strengths and
  weaknesses: check the results.
