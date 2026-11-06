"""The record and the summary of the gaps of the collection of recipes
(phase 6-ter, PR 5b, decision 6)."""
import json
import logging
import sys

import pytest

from astrocook import __version__, gaps


def lines(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def test_the_file_by_system(monkeypatch, tmp_path):
    monkeypatch.delenv(gaps.ENV_VAR)
    monkeypatch.setattr("pathlib.Path.home", classmethod(lambda cls: tmp_path))
    monkeypatch.delenv("XDG_DATA_HOME", raising=False)
    monkeypatch.setattr(sys, "platform", "darwin")
    assert gaps.gaps_path() == (tmp_path / "Library" / "Application Support"
                                / "Astrocook" / "gaps.jsonl")
    monkeypatch.setattr(sys, "platform", "linux")
    assert gaps.gaps_path() == tmp_path / ".local" / "share" / "astrocook" / "gaps.jsonl"
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    assert gaps.gaps_path() == tmp_path / "xdg" / "astrocook" / "gaps.jsonl"
    monkeypatch.setenv("XDG_DATA_HOME", "relative/dir")  # not absolute: ignored
    assert gaps.gaps_path() == tmp_path / ".local" / "share" / "astrocook" / "gaps.jsonl"
    monkeypatch.setattr(sys, "platform", "win32")
    monkeypatch.setenv("APPDATA", str(tmp_path / "roaming"))
    assert gaps.gaps_path() == tmp_path / "roaming" / "Astrocook" / "gaps.jsonl"
    monkeypatch.delenv("APPDATA")
    assert gaps.gaps_path() == (tmp_path / "AppData" / "Roaming" / "Astrocook"
                                / "gaps.jsonl")


def test_the_variable_moves_or_turns_off_the_record(monkeypatch, tmp_path):
    monkeypatch.setenv(gaps.ENV_VAR, str(tmp_path / "other" / "g.jsonl"))
    assert gaps.gaps_path() == tmp_path / "other" / "g.jsonl"
    assert gaps.record_gap("declared", "x", "s")
    assert (tmp_path / "other" / "g.jsonl").is_file()  # the folder is made
    monkeypatch.setenv(gaps.ENV_VAR, "OFF")
    assert gaps.gaps_path() is None
    assert not gaps.record_gap("declared", "x", "s")
    assert gaps.read_gaps() == ([], 0)


def test_a_line_holds_only_five_fields(tmp_path):
    path = gaps.gaps_path()
    assert gaps.record_gap("declared", "Needed  z_em.\nMeasured it.", "mock")
    assert gaps.record_gap("failed", "x" * 5000, "mock")
    first, second = lines(path)
    assert set(first) == {"date", "version", "kind", "text", "session"}
    assert first["version"] == __version__ and first["session"] == "mock"
    assert first["kind"] == "declared" and first["text"] == "Needed z_em. Measured it."
    assert first["date"].endswith("+00:00")
    assert len(second["text"]) == gaps.MAX_TEXT
    with pytest.raises(ValueError, match="Unknown kind"):
        gaps.record_gap("other", "x", "s")


def test_a_record_that_cannot_be_written_does_not_raise(monkeypatch, tmp_path, caplog):
    blocked = tmp_path / "file"
    blocked.write_text("")
    monkeypatch.setenv(gaps.ENV_VAR, str(blocked / "sub" / "g.jsonl"))  # parent is a file
    monkeypatch.setattr(gaps, "_warned", False)
    with caplog.at_level(logging.WARNING, logger="astrocook.gaps"):
        assert not gaps.record_gap("declared", "x", "s")
        assert not gaps.record_gap("declared", "y", "s")
    assert len([r for r in caplog.records if "cannot be recorded" in r.message]) == 1


def test_read_skips_the_lines_it_cannot_read():
    path = gaps.gaps_path()
    assert gaps.read_gaps() == ([], 0)  # no file: no gaps
    gaps.record_gap("declared", "one", "s")
    with open(path, "a") as f:
        f.write("not json\n\n[1]\n" + json.dumps({"date": "2026-01-01"}) + "\n")
    gaps.record_gap("failed", "two", "s")
    records, skipped = gaps.read_gaps()
    assert [r["text"] for r in records] == ["one", "two"] and skipped == 3


def rec(kind, text, session="s", date="2026-11-01T10:00:00+00:00", version="2.5.0"):
    return dict(date=date, version=version, kind=kind, text=text, session=session)


def test_summary_groups_by_kind_and_text_without_numbers():
    records = [
        rec("prerequisite", "flux_x needs z_em."),
        rec("prerequisite", "flux_x needs z_em.", session="t", date="2026-11-05T10:00:00+00:00",
            version="2.6.0"),
        rec("declared", "Needed z_em: Ly-a peak at 2.345."),
        rec("declared", "needed  z_em: ly-a peak at 3.1.", session="t"),
        rec("declared", "Needed a mask."),
        rec("undone", "absorbers_doublets_auto undone."),
    ]
    groups = gaps.summarize(records)
    assert [(g["kind"], g["count"]) for g in groups] == [
        ("declared", 2), ("prerequisite", 2), ("declared", 1), ("undone", 1)]
    g = groups[0]
    assert g["text"] == "Needed z_em: Ly-a peak at 2.345." and g["sessions"] == 2
    assert g["examples"] == ["Needed z_em: Ly-a peak at 2.345.",
                             "needed  z_em: ly-a peak at 3.1."]
    pre = groups[1]
    assert (pre["first"], pre["last"], pre["versions"], pre["sessions"]) == (
        "2026-11-01T10:00:00+00:00", "2026-11-05T10:00:00+00:00", ["2.5.0", "2.6.0"], 2)
    assert [g["kind"] for g in gaps.summarize(records, kind="undone")] == ["undone"]
    assert sum(g["count"] for g in gaps.summarize(records, since="2026-11-02")) == 1
    assert len(gaps.summarize(records, limit=2)) == 2
    with pytest.raises(ValueError, match="not a date"):
        gaps.summarize(records, since="yesterday")
    with pytest.raises(ValueError, match="Unknown kind"):
        gaps.summarize(records, kind="other")


def test_the_text_of_the_summary():
    assert gaps.format_summary([]) == "No gaps recorded."
    text = gaps.format_summary(gaps.summarize([
        rec("declared", "Needed A."), rec("declared", "Needed A.", session="t"),
        rec("undone", "x undone.")]), skipped=2)
    assert text.splitlines()[0] == "declared (2):"
    assert "2x in 2 session(s), 2026-11-01 to 2026-11-01, versions 2.5.0" in text
    assert "undone (1):" in text and "(2 line(s)" in text


def test_the_command(capsys, tmp_path):
    assert gaps.main(["--path"]) == 0
    assert capsys.readouterr().out.strip() == str(gaps.gaps_path())
    assert gaps.main([]) == 0
    assert "No gaps recorded." in capsys.readouterr().out
    gaps.record_gap("declared", "Needed z_em.", "mock")
    gaps.record_gap("failed", "t failed on a valid request: TypeError.", "mock")
    assert gaps.main(["--kind", "failed"]) == 0
    out = capsys.readouterr().out
    assert "failed (1):" in out and "declared" not in out
    assert gaps.main(["--json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["skipped_lines"] == 0 and len(data["groups"]) == 2
    with pytest.raises(SystemExit):
        gaps.main(["--since", "yesterday"])
    other = tmp_path / "other.jsonl"
    assert gaps.main(["--file", str(other)]) == 0
    assert "No gaps recorded." in capsys.readouterr().out


def test_the_command_with_the_record_off(monkeypatch, capsys):
    monkeypatch.setenv(gaps.ENV_VAR, "off")
    assert gaps.main([]) == 0
    assert "is off" in capsys.readouterr().out
    assert gaps.main(["--path"]) == 0
    assert capsys.readouterr().out.startswith("off")
