"""The reader of the progress lines of the recipes (phase 6-ter, PR 4)."""
import logging
import threading

import pytest

from astrocook.recipes.progress import (
    PROGRESS_MARKER, is_progress_line, parse_progress, watch_progress)

logger = logging.getLogger("astrocook.test_recipes_progress")


@pytest.mark.parametrize("message, value", [
    (">> PROGRESS: 42", 42),
    ("INFO >> PROGRESS: 7\n", 7),      # as formatted by the GUI handler
    (">> PROGRESS: many", None),
    (">> PROGRESS:", None),
    ("Fitting group 3 of 7", None)])
def test_parse_progress(message, value):
    assert parse_progress(message) == value


def test_a_progress_line_with_an_unreadable_value_is_still_a_progress_line():
    assert is_progress_line(f"{PROGRESS_MARKER} many")
    assert not is_progress_line("Fitting group 3 of 7")


def test_watch_progress_gives_the_percent_and_the_last_status():
    seen = []
    with watch_progress(lambda p, s: seen.append((p, s))):
        logger.info("Starting")
        logger.info(">> PROGRESS: 0")
        logger.info("Group 1")
        logger.info("Group 2")
        logger.info(">> PROGRESS: 50")
        logger.info(">> PROGRESS: x")  # unreadable: ignored
    assert seen == [(0, "Starting"), (50, "Group 2")]


def test_watch_progress_lowers_the_level_and_restores_it():
    package = logging.getLogger("astrocook")
    level, handlers = package.level, list(package.handlers)
    package.setLevel(logging.WARNING)
    try:
        # The root logger of the application is at WARNING: the lines reach
        # the handler only while the block runs
        seen = []
        quiet = logging.getLogger("astrocook.test_recipes_progress.quiet")
        quiet.setLevel(logging.NOTSET)
        with watch_progress(lambda p, s: seen.append(p)):
            quiet.info(">> PROGRESS: 10")
        quiet.info(">> PROGRESS: 20")
        assert seen == [10]
        assert package.level == logging.WARNING
        assert package.handlers == handlers
    finally:
        package.setLevel(level)


def test_watch_progress_restores_the_logger_after_an_error():
    package = logging.getLogger("astrocook")
    handlers = list(package.handlers)
    with pytest.raises(RuntimeError), watch_progress(lambda p, s: None):
        raise RuntimeError("boom")
    assert package.handlers == handlers


def test_watch_progress_ignores_other_threads():
    seen = []
    with watch_progress(lambda p, s: seen.append(p)):
        t = threading.Thread(target=lambda: logger.info(">> PROGRESS: 99"))
        t.start()
        t.join()
        logger.info(">> PROGRESS: 1")
    assert seen == [1]
