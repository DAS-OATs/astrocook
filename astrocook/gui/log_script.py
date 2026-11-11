"""The session log as the text of the Log Scripter, and back.

Each entry is a call, ``recipe(param=value, ...)``; the entries undone are
commented out. The reason of an entry is a comment before it
(``# reason: ...``), a note a comment where it was written
(``# note: ...``); the lines after the first of a long text start with
``#`` and spaces. *Run All* skips the comments as steps, but keeps the
reasons and the notes in the log it rebuilds (phase 6-ter, PR 5a): a
replay does not lose them.
"""
import re
from typing import List, Optional, Tuple

from astrocook.core.structures import NOTE_KINDS, HistoryLog, LogEntry

#: The tags of the comments that hold a reason or a note
TAGS = ("reason", *NOTE_KINDS)
_TAG_LINE = re.compile(r"#\s?(" + "|".join(TAGS) + r")\s*:\s?(.*)$")
_MORE_LINE = re.compile(r"#\s{2,}(.*)$")


def _tagged(tag: str, text: str) -> List[str]:
    head = f"# {tag}: "
    lines = text.splitlines() or [""]
    return [head + lines[0]] + ["#" + " " * (len(head) - 1) + line for line in lines[1:]]


def _call(entry: LogEntry) -> str:
    params = ", ".join(f"{k}='{v}'" if isinstance(v, str) else f"{k}={v}"
                       for k, v in entry.params.items())
    return f"{entry.recipe_name}({params})"


#: The first line of the log of the archive a session was opened from, in
#: the Log Scripter: its lines start with ``##``, which Run All skips and
#: does not read as reasons or notes (F109)
PREVIOUS_HEADER = ("## The log of the archive this session was opened from: read only, "
                   "not run by Run All, not saved again.")


def log_to_script(log: HistoryLog) -> List[str]:
    """The lines of the Log Scripter for ``log``, the entries undone and
    their notes commented out; first, the log of the archive the session
    was opened from, if any, as ``##`` comments."""
    out = []
    previous = getattr(log, 'previous', None)
    if previous is not None:
        out += [PREVIOUS_HEADER] + [f"## {line}" for line in log_to_script(previous)] + [""]
    for item, active in log.timeline(include_undone=True):
        if isinstance(item, LogEntry):
            lines = (_tagged("reason", item.reason) if item.reason else []) + [_call(item)]
        else:
            lines = _tagged(item.kind, item.text)
        out += lines if active else [f"# {line}" for line in lines]
    return out


class ScriptComments:
    """Reads the reasons and the notes of a script while it runs: each
    comment line goes to :meth:`read`, and before each call :meth:`notes`
    and :meth:`reason` give what was written before it."""

    def __init__(self):
        self._tag: Optional[str] = None
        self._lines: List[str] = []
        self._notes: List[Tuple[str, str]] = []
        self._reason: List[str] = []

    def _close(self):
        if self._tag == "reason":
            self._reason = self._lines
        elif self._tag is not None:
            self._notes.append((self._tag, "\n".join(self._lines)))
        self._tag, self._lines = None, []

    def read(self, line: str) -> bool:
        """Read a line of the script (stripped); True if it was a reason, a
        note or one of their lines."""
        tagged = _TAG_LINE.match(line)
        if tagged:
            self._close()
            self._tag, self._lines = tagged.group(1), [tagged.group(2)]
            return True
        more = _MORE_LINE.match(line)
        if more and self._tag is not None:
            self._lines.append(more.group(1))
            return True
        self._close()
        return False

    def notes(self) -> List[Tuple[str, str]]:
        """The notes read since the last call, as ``(kind, text)``."""
        self._close()
        notes, self._notes = self._notes, []
        return [(kind, text) for kind, text in notes if text.strip()]

    def reason(self) -> Optional[str]:
        """The reason read for the next call, if any."""
        self._close()
        reason, self._reason = self._reason, []
        return "\n".join(reason) or None
