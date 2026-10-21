# Using Astrocook from an AI Agent

:::{warning}
This is a **pilot**: a few recipes only, to learn what an agent needs
from Astrocook. Tools, parameters and outputs may change in later
versions.
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

Each tool is named after the recipe it runs, prefixed with its menu, and
uses the same defaults as the recipe dialog of the GUI.

| Tool | What it does |
|---|---|
| `session_open` | Opens a spectrum or a session file (`.acs2`, V1 `.acs`) as a new session; optionally sets z_em and the resolving power (*Edit > Set Properties*). |
| `continuum_estimate_auto` | *Continuum > Auto-estimate Continuum*. Returns the fraction of absorbed pixels (in the Lyman-alpha forest and on the red side) and the S/N. |
| `absorbers_identify_lines` | *Absorbers > Identify Absorption Lines*, without creating components. Returns the identified regions with their candidate transitions and scores, and the candidate systems. |
| `absorbers_doublets_auto` | *Absorbers > Auto-detect Doublets*: identifies and fits the CIV, SiIV and MgII doublets on the red side. Returns the new components, with errors, reduced chi2 and flags for the ones stuck at a bound of the fit. |
| `session_save` | Saves the session as `.acs2`, with the original data and the session log. |

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

If something goes wrong (a missing file, z_em not set), the tool returns
Astrocook's error message, and the agent can usually correct its request
by itself.

## 5. Reviewing the Result in the GUI

Open the saved archive in the GUI (`astrocook /data/J0103-1305_agent.acs2`).
The session log (*View > View Session Log*) lists the recipes the agent ran, with
their parameters, as if you had run them from the menus: you can inspect
the components in the System Inspector, correct them, refit them, or
replay the log as a script.

## 6. Limits of the Pilot

- Only the five tools above: no manual editing of components, no Lyman
  forest fit, no export.
- Paths are on the computer where the server runs; give them in full.
- `absorbers_doublets_auto` removes the components already on the red
  side before fitting, and replaces the absorbed mask with one made with
  the default parameters. A call can take minutes on a large spectrum.
- The recipes are the same as in the GUI, with the same strengths and
  weaknesses: check the results.
