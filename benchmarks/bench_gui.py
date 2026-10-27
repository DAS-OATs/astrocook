"""Timing benchmark of the GUI, to compare before and after a change.

Not a test (not collected by pytest, not run in CI): run it by hand on the
same machine before and after a change that touches the GUI, and compare the
numbers. The GUI runs offscreen on the small J0103-1305 archive of the tests
(``tests/data/J0103-1305_metallines_lya_new_3.acs2``: 2156 pixels, 764
components, no continuum). Five cases:

- ``inspector_first``: the first opening of the System Inspector in a new
  main window (the inspector is built, the table filled, the first row
  selected, the main plot highlights the selection and the velocity plot
  draws it; median of ``--repeat``, a new main window each time);
- ``inspector_again``: hiding it and opening it again (median of
  ``--repeat``);
- ``main_plot``: a full redraw of the main plot (median of ``--repeat``);
- ``refit``: *Refit Selected* on the first CIV component, from the
  request to the end of the result processing (recipe in the worker thread,
  new state in the history, plots and table updated; median of
  ``--repeat``). The archive has no continuum: a first refit, not timed,
  estimates it, so that the timed ones do not. Split into ``refit_recipe``
  (the recipe in the worker thread) and ``refit_gui`` (the rest: the
  difference of the medians);
- ``select_row``: selecting another row of the inspector table (the main
  plot highlights the new selection, the velocity plot shows it; median of
  ``--repeat``, a different row each time). Run last, with the main window
  at 1440 x 860 and the inspector at 1000 x 700 (a laptop screen; the cost
  of a draw grows with the size of the plots), after the refits (the
  session has a continuum and a model).

Every case ends after the pending Qt events are processed, so that the
draws deferred with ``draw_idle`` are included.

Usage::

    python benchmarks/bench_gui.py [--repeat N] [--json out.json]
                                   [--profile DIR]

``--profile DIR`` also writes a ``cProfile`` file per case into ``DIR``
(``<case>.prof``, to read with ``pstats`` or ``snakeviz``); the profiler
slows the run down, so the times of a profiled run are not comparable with
the others.

The GUI runs offscreen unless ``QT_QPA_PLATFORM`` is set: e.g.
``QT_QPA_PLATFORM=cocoa`` on macOS shows the windows, and the times then
include the drawing on the real screen.

A run takes about 2.5 minutes with the default ``--repeat`` (5). The times
depend on the machine: only compare runs made on the same one.
"""
import argparse
import cProfile
import json
import logging
import os
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SESSION_FILE = ROOT / "tests" / "data" / "J0103-1305_metallines_lya_new_3.acs2"


def _setup_qt():
    """Offscreen QApplication, settings kept away from the user's ones, and
    message boxes that do not block."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ["XDG_CONFIG_HOME"] = tempfile.mkdtemp(prefix="bench_gui_")
    from PySide6.QtWidgets import QApplication, QMessageBox
    app = QApplication.instance() or QApplication([])
    for kind in ("warning", "critical", "information", "question"):
        setattr(QMessageBox, kind,
                staticmethod(lambda *a, **k: QMessageBox.StandardButton.Ok))
    return app


def _flush(app):
    """Process the pending events (deferred draws, queued signals)."""
    for _ in range(3):
        app.processEvents()


class Bench:
    def __init__(self, app, profile_dir=None):
        self.app = app
        self.profile_dir = profile_dir
        self.window = None
        self.worker_times = []
        self._time_the_workers()

    def _time_the_workers(self):
        """Record how long each recipe runs in its worker thread, to separate
        the time of the recipe from that of the GUI."""
        from astrocook.gui.qt_workers import RecipeWorker
        run = RecipeWorker.run

        def timed_run(worker):
            t0 = time.perf_counter()
            try:
                run(worker)
            finally:
                self.worker_times.append(time.perf_counter() - t0)

        RecipeWorker.run = timed_run

    def new_window(self):
        # Imported after the QApplication exists: the GUI selects the QtAgg
        # matplotlib backend at import time (F7)
        from astrocook.gui.main_window import MainWindowV2
        if self.window is not None:
            self.window.close()
            self.window.deleteLater()
            _flush(self.app)
        win = MainWindowV2(None, None)
        win._show_custom_message = lambda *a, **k: None
        win.show()
        win.open_session_from_path(str(SESSION_FILE))
        _flush(self.app)
        self.window = win
        return win

    def time(self, name, func, repeat, setup=None):
        """Median wall time of ``func`` over ``repeat`` runs, each after
        ``setup`` (not timed); profiled into ``<name>.prof`` if asked."""
        profiler = cProfile.Profile() if self.profile_dir else None
        times = []
        for _ in range(repeat):
            if setup:
                setup()
                _flush(self.app)
            if profiler:
                profiler.enable()
            t0 = time.perf_counter()
            func()
            _flush(self.app)
            times.append(time.perf_counter() - t0)
            if profiler:
                profiler.disable()
        if profiler:
            profiler.dump_stats(str(Path(self.profile_dir) / f"{name}.prof"))
        return float(np.median(times))

    def open_inspector(self):
        self.window._on_view_system_inspector()

    def reopen_inspector(self):
        self.window.system_inspector.hide()
        _flush(self.app)
        self.window._on_view_system_inspector()

    def redraw_main_plot(self):
        viewer = self.window.plot_viewer
        viewer.plot_spectrum(self.window.active_history.current_state,
                             initial_draw=True)

    def select_civ(self):
        """Select the first CIV component of the archive (z = 3.0085), a
        small group outside the forest: the refit is short and its time is
        mostly that of the GUI."""
        systs = self.window.active_history.current_state.systs
        comp = next(c for c in systs.components if c.series == "CIV")
        self.window.system_inspector.focus_on_component(comp.uuid)

    def laptop_size(self):
        self.window.resize(1440, 860)
        self.window.system_inspector.resize(1000, 700)
        self._row = 10

    def select_next_row(self):
        self._row += 7
        self.window.system_inspector.table_view.selectRow(self._row)

    def refit_selected(self):
        inspector = self.window.system_inspector
        comps = inspector._get_selected_components()
        assert comps, "no component selected"
        before = self.window.active_history.current_state
        inspector._refit_list(comps)
        # The result is processed after the worker ends (with a QTimer):
        # wait for the new state in the history
        while self.window.active_history.current_state is before:
            self.app.processEvents()
            time.sleep(0.005)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--repeat", type=int, default=5,
                        help="repetitions of the repeated cases")
    parser.add_argument("--json", help="write the results to this file")
    parser.add_argument("--profile", metavar="DIR",
                        help="write a cProfile file per case into DIR")
    args = parser.parse_args()
    logging.basicConfig(level=logging.WARNING)
    if args.profile:
        Path(args.profile).mkdir(parents=True, exist_ok=True)

    app = _setup_qt()
    bench = Bench(app, args.profile)
    results = {"inspector_first_s": bench.time(
        "inspector_first", bench.open_inspector, args.repeat,
        setup=bench.new_window)}
    results["inspector_again_s"] = bench.time(
        "inspector_again", bench.reopen_inspector, args.repeat)
    results["main_plot_s"] = bench.time(
        "main_plot", bench.redraw_main_plot, args.repeat)
    bench.select_civ()
    _flush(app)
    bench.refit_selected()  # estimates the continuum, not timed
    _flush(app)
    bench.worker_times.clear()
    results["refit_s"] = bench.time("refit", bench.refit_selected,
                                    args.repeat)
    results["refit_recipe_s"] = float(np.median(bench.worker_times))
    results["refit_gui_s"] = results["refit_s"] - results["refit_recipe_s"]
    bench.laptop_size()
    _flush(app)
    results["select_row_s"] = bench.time("select_row", bench.select_next_row,
                                         args.repeat)
    for key, value in results.items():
        print(f"{key:20s} {value:10.3f}")
    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
