# Prompts of the trials

Given to each agent as written (with `<folder>` replaced), in a fresh context:
no knowledge of Astrocook's code, tests or the true parameters of the mocks.
Phase 6-ter, PR 7 (`ROADMAP.md`, section 4.3).

## Rules (every trial)

> You are helping an astronomer analyse a quasar spectrum with Astrocook, a
> spectral-analysis package that you can use ONLY through its MCP server. You
> know nothing about Astrocook's code, and you must keep it that way: this is a
> usability trial of the server, and it is only valid if you see exactly what an
> MCP client shows a model.
>
> - Work only in the folder `<folder>` (cd there first).
> - Do NOT read, list, grep or open anything under the repository, any Python
>   package source (site-packages included), any file outside your folder, or
>   the client script that the `./mcp` wrapper calls. Do not run `docker`
>   yourself, do not inspect processes. Do not search the web for Astrocook. If
>   you are unsure whether something is allowed, don't.
> - You talk to the server only with the wrapper in your folder: `./mcp tools`
>   prints the server's instructions and every tool (name, description, JSON
>   input schema), exactly as a client gives them to a model (run it first and
>   read it); `./mcp call <tool_name> '<JSON arguments>'` calls a tool and prints
>   its result. Images are saved as PNG files and the result gives their path:
>   open them with your file-reading tool to see them.
> - (Docker) The server runs in a container: it sees your folder `./data` as
>   `/data`. (Local) The server runs on this machine: give it absolute paths.
> - Some tools may run for minutes: that is expected, wait for them.
> - You may use your own Python/shell on files in your folder if the server
>   cannot do something, but every time you work outside the server's tools,
>   declare it as a gap with the tool `session_note` and `kind="gap"` (what you
>   needed, what you did instead). Also declare as a gap anything you needed
>   that no tool provided.
> - Budget: at most about 60 tool calls of the server in total. Be economical
>   with large outputs.
>
> When you are done, write `report.md` in your folder with two parts, and also
> return its full text as your final answer: (1) **Answer** to the astronomer;
> (2) **Critique of the server**, honest and concrete: what worked; what was
> confusing or wrong (tool descriptions, parameters, outputs, error messages);
> what was missing (gaps); how easy it was to find the right tool among all
> those offered (did the number of tools get in the way? would fewer, grouped
> tools be better or worse?); for the slow calls, how long they took and whether
> you would have wanted progress information, an estimate of the time, or a way
> to cancel them; anything you did outside the server and why.

## Trial 1: the full-chain mock (`make_data.py mock`, Docker)

> "I have a simulated quasar spectrum, /data/spectrum.acs2. Please find the
> metal absorption systems on the red side of the Lyman-alpha emission, fit
> them, and tell me which components are reliable and which are not, and why.
> I'd like to be able to reproduce what you did, so leave me a saved session
> and a description of the methods."

## Trial 2: J0103-1305 (`make_data.py j0103`, local)

> "This is a UVES spectrum of the quasar J0103-1305, as it comes from the
> reduction: nothing has been done on it yet. Please estimate the continuum,
> tell me which metal absorption systems there are, fit them, and tell me which
> results I can trust and which I can't, and why. I don't want to lose any work
> along the way, and I need to be able to reproduce what you did: leave me a
> saved session, a log I could replay, and a description of the methods for a
> paper."

## Trial 3: editing the components (`make_data.py edit`, local)

> "In this simulated spectrum, /data/spectrum.acs2, the CIV system near
> z = 2.20 is the one I care about. Fit it, then check whether the model is
> adequate: try the alternatives you think are sensible (more or fewer
> components, constraints between the ions of the system), keep the model you
> can defend and explain why. If a change makes things worse, go back to the
> model before. At the end tell me the final components, what you tried and
> rejected, and leave me a saved session and a description of the methods."

Run locally (PR 7, after PR 7-bis), with the path of the file in the rules
and one line after the request: the server runs locally, so
`/data/spectrum.acs2` is the file in the folder's `data/`.
