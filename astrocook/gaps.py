"""The gaps of the collection of recipes (phase 6-ter, PR 5b, decision 6).

A *gap* is something the analysis needed and no recipe provided. The gaps
are recorded where they are met (by the MCP server today, by any caller
that wants to) in a local file of the user, one JSON object per line, and
summarized by :func:`summarize` (the command ``astrocook-gaps``, the tool
``gaps_summary`` of the server) as the input of the plan of a phase.

Privacy: the file is local and **never sent anywhere**. A line holds only
``date``, ``version``, ``kind``, ``text`` and ``session`` (the name of the
session), nothing else of the user: no paths, no parameter values, no data
of the spectrum. The text is what the agent wrote (kind ``declared``) or a
fixed sentence of the server naming a tool, a prerequisite or an exception
type (never its message, which may hold a path).

Where the file is (``ASTROCOOK_GAPS_FILE`` overrides it; ``off`` disables
the record):

=========  =======================================================
macOS      ``~/Library/Application Support/Astrocook/gaps.jsonl``
Linux      ``$XDG_DATA_HOME/astrocook/gaps.jsonl``, else
           ``~/.local/share/astrocook/gaps.jsonl``
Windows    ``%APPDATA%\\Astrocook\\gaps.jsonl``
=========  =======================================================
"""
import argparse
import json
import logging
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

logger = logging.getLogger(__name__)

#: The kinds of gap. ``declared``: written by the agent (what was needed and
#: what it did instead); ``prerequisite``: a prerequisite that no recipe
#: estimates from the data; ``failed``: a recipe that failed with an
#: unexpected error on a valid request; ``undone``: the result of a recipe
#: that the agent undid.
KINDS = ("declared", "prerequisite", "failed", "undone")

#: The environment variable that moves the file, or turns the record off
ENV_VAR = "ASTROCOOK_GAPS_FILE"

#: The longest text recorded; the rest is cut
MAX_TEXT = 2000

_warned = False


def gaps_path() -> Optional[Path]:
    """The file of the gaps, or None if the record is off."""
    given = os.environ.get(ENV_VAR)
    if given is not None and given.strip():
        if given.strip().lower() == "off":
            return None
        return Path(given).expanduser()
    home = Path.home()
    if sys.platform == "darwin":
        return home / "Library" / "Application Support" / "Astrocook" / "gaps.jsonl"
    if sys.platform.startswith("win"):
        base = os.environ.get("APPDATA") or str(home / "AppData" / "Roaming")
        return Path(base) / "Astrocook" / "gaps.jsonl"
    data = os.environ.get("XDG_DATA_HOME")
    base = Path(data) if data and os.path.isabs(data) else home / ".local" / "share"
    return base / "astrocook" / "gaps.jsonl"


def record_gap(kind: str, text: str, session: str) -> bool:
    """Append a gap to the file of the user.

    Never raises: a record that cannot be written must not break the caller
    (one warning is logged). Returns True if the line was written.

    Parameters
    ----------
    kind : str
        One of :data:`KINDS`.
    text : str
        What was needed and what was done instead, or the sentence of the
        server; cut to :data:`MAX_TEXT` characters.
    session : str
        The name of the session.
    """
    global _warned
    if kind not in KINDS:
        raise ValueError(f"Unknown kind of gap '{kind}': one of {', '.join(KINDS)}.")
    path = gaps_path()
    if path is None:
        return False
    from astrocook import __version__
    line = json.dumps({
        "date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "version": __version__, "kind": kind,
        "text": " ".join(str(text).split())[:MAX_TEXT], "session": str(session)},
        ensure_ascii=False)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "a", encoding="utf-8") as f:
            f.write(line + "\n")
    except OSError as e:
        if not _warned:
            _warned = True
            logger.warning(f"The gaps cannot be recorded in {path}: {e}")
        return False
    return True


def read_gaps(path: Optional[Path] = None) -> Tuple[List[Dict[str, Any]], int]:
    """The gaps recorded, and the number of lines skipped (not JSON, or
    without the fields). A file that does not exist holds no gaps."""
    path = path if path is not None else gaps_path()
    if path is None or not Path(path).is_file():
        return [], 0
    records, skipped = [], 0
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            try:
                r = json.loads(line)
                ok = (isinstance(r, dict)
                      and all(isinstance(r.get(k), str) for k in
                              ("date", "version", "kind", "text", "session")))
            except ValueError:
                ok = False
            if ok:
                records.append(r)
            else:
                skipped += 1
    return records, skipped


def _normal(text: str) -> str:
    """The text in lower case, with the spaces reduced and the numbers
    masked: the key that groups the gaps."""
    return re.sub(r"\d+(?:\.\d+)?", "#", " ".join(text.lower().split()))


def summarize(records: Iterable[Dict[str, Any]], since: Optional[str] = None,
              kind: Optional[str] = None, limit: Optional[int] = None
              ) -> List[Dict[str, Any]]:
    """Group the gaps by kind and text (without the numbers).

    The texts of the kind ``declared`` are free, so their groups are only
    those of equal wording: the reader joins the similar ones.

    Parameters
    ----------
    records : iterable of dict
        The output of :func:`read_gaps`.
    since : str, optional
        Only the gaps from this date on (``YYYY-MM-DD``).
    kind : str, optional
        Only this kind.
    limit : int, optional
        At most this number of groups, the most frequent first.

    Returns
    -------
    list of dict
        Each with ``kind``, ``text`` (the first one met), ``count``,
        ``sessions`` (how many different), ``versions``, ``first``,
        ``last`` (dates) and ``examples`` (up to three different texts),
        the most frequent first.
    """
    if kind is not None and kind not in KINDS:
        raise ValueError(f"Unknown kind of gap '{kind}': one of {', '.join(KINDS)}.")
    if since is not None:
        try:
            datetime.strptime(since, "%Y-%m-%d")
        except ValueError:
            raise ValueError(f"since: '{since}' is not a date (YYYY-MM-DD).") from None
    groups: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for r in records:
        if kind is not None and r["kind"] != kind:
            continue
        if since is not None and r["date"][:10] < since:
            continue
        g = groups.setdefault((r["kind"], _normal(r["text"])), {
            "kind": r["kind"], "text": r["text"], "count": 0, "_sessions": set(),
            "_versions": set(), "first": r["date"], "last": r["date"], "examples": []})
        g["count"] += 1
        g["_sessions"].add(r["session"])
        g["_versions"].add(r["version"])
        g["first"], g["last"] = min(g["first"], r["date"]), max(g["last"], r["date"])
        if r["text"] not in g["examples"] and len(g["examples"]) < 3:
            g["examples"].append(r["text"])
    order = {k: i for i, k in enumerate(KINDS)}
    out = []
    for g in groups.values():
        g["sessions"] = len(g.pop("_sessions"))
        g["versions"] = sorted(g.pop("_versions"))
        out.append(g)
    out.sort(key=lambda g: (-g["count"], order.get(g["kind"], len(order)), g["text"]))
    return out[:limit] if limit is not None else out


def format_summary(groups: List[Dict[str, Any]], skipped: int = 0) -> str:
    """The groups of :func:`summarize` as text, by kind."""
    if not groups:
        return "No gaps recorded."
    lines = []
    for kind in KINDS:
        of_kind = [g for g in groups if g["kind"] == kind]
        if not of_kind:
            continue
        lines.append(f"{kind} ({sum(g['count'] for g in of_kind)}):")
        for g in of_kind:
            lines.append(f"  {g['count']}x in {g['sessions']} session(s), "
                         f"{g['first'][:10]} to {g['last'][:10]}, "
                         f"versions {', '.join(g['versions'])}")
            lines.append(f"    {g['text']}")
            for example in g["examples"][1:]:
                lines.append(f"    also: {example}")
    if skipped:
        lines.append(f"({skipped} line(s) of the file could not be read)")
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    """Entry point of ``astrocook-gaps``: summarize the gaps recorded."""
    parser = argparse.ArgumentParser(
        prog="astrocook-gaps",
        description="Summarize the gaps of the collection of recipes: what the "
                    "analyses needed and no recipe provided, recorded in a local "
                    f"file (never sent anywhere; {ENV_VAR}=off disables it).")
    parser.add_argument("--since", metavar="YYYY-MM-DD", help="only the gaps from this date")
    parser.add_argument("--kind", choices=KINDS, help="only this kind")
    parser.add_argument("--json", action="store_true", help="print the groups as JSON")
    parser.add_argument("--path", action="store_true",
                        help="print the path of the file and exit")
    parser.add_argument("--file", type=Path, help="read this file instead")
    args = parser.parse_args(argv)
    path = args.file if args.file is not None else gaps_path()
    if args.path:
        print(path if path is not None else f"off ({ENV_VAR}=off)")
        return 0
    if path is None:
        print(f"The record of the gaps is off ({ENV_VAR}=off).")
        return 0
    records, skipped = read_gaps(path)
    try:
        groups = summarize(records, since=args.since, kind=args.kind)
    except ValueError as e:
        parser.error(str(e))
    if args.json:
        print(json.dumps({"path": str(path), "skipped_lines": skipped, "groups": groups},
                         indent=2, ensure_ascii=False))
    else:
        print(f"{path}\n{format_summary(groups, skipped)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
