"""Rebuild the screenshots of the tutorials from the GUI.

Each screenshot in ``docs/_static`` is a *scene*: the spectra to open, the
recipes to run, the state of the view (axis limits, panels, selections, an
open menu or dialog) and the part of the window to keep. The scenes run the
real GUI offscreen, at twice the logical resolution (like a retina screen),
so the images can be rebuilt after every change of the interface:

    python docs/screenshots/make_screenshots.py            # all scenes
    python docs/screenshots/make_screenshots.py --only continuum_auto
    python docs/screenshots/make_screenshots.py --list

The spectra ("My QSO" and "My other QSO") are read from ``--data``
(default ``docs/screenshots/data``, not in the repository; ``--download``
fetches them); see ``docs/screenshots/README.md``. Windows get a neutral frame (title bar,
rounded corners, shadow) and the mouse pointer is drawn where a scene
puts it: the offscreen platform has neither.

Look at every image before committing it: the scenes check that they run,
not that the result reads well.
"""
import argparse
import io
import json
import os
import shutil
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

# The platform, the scale and the settings must be fixed before Qt starts.
# A fresh settings directory gives every run the default state of the GUI
# (sections of the plot controls, last folders).
# The offscreen platform has a screen of 400 x 400 points by default, which
# squeezes the dialogs: a large screen, scaled by two (like a retina one).
_settings = tempfile.mkdtemp(prefix="astrocook-shots-")
_screen = Path(_settings) / "screen.json"
_screen.write_text(json.dumps({"screens": [{
    "name": "docs", "x": 0, "y": 0, "width": 6400, "height": 4000,
    "logicalDpi": 96, "logicalBaseDpi": 96, "dpr": 1}]}))
os.environ["QT_QPA_PLATFORM"] = f"offscreen:configfile={_screen}"
os.environ["QT_SCALE_FACTOR"] = "2"
os.environ["XDG_CONFIG_HOME"] = _settings

from PIL import Image, ImageChops, ImageDraw, ImageFilter  # noqa: E402
from PySide6.QtCore import QBuffer, QIODevice, QPoint, QRect, Qt, QTimer  # noqa: E402
from PySide6.QtGui import QColor, QFont, QImage, QPainter  # noqa: E402
from PySide6.QtWidgets import QApplication, QMenu, QPushButton, QToolTip  # noqa: E402

SCALE = 2                       # physical pixels per logical pixel
WINDOW_SIZE = (1000, 640)       # logical size of the main window
TITLE_BAR = 28                  # logical height of the drawn title bar
RADIUS = 10                     # logical radius of the window corners
SHADOW = 24                     # logical margin left for the shadow
HERE = Path(__file__).resolve().parent
DEFAULT_DATA = HERE / "data"
DEFAULT_OUT = HERE.parent / "_static"
SESSIONS = ("My QSO", "My other QSO")
# The spectra are assets of a pre-release, downloaded on the first run.
# The session takes the name of the file: underscores become spaces.
DATA_RELEASE = "https://github.com/DAS-OATs/astrocook/releases/download/data-docs-2026"
# Properties of "My QSO", as set in the tutorial on session properties
Z_EM = 2.730
RESOL = 135000

SCENES = {}


def scene(func):
    """Register a scene; its name is the name of the image."""
    SCENES[func.__name__] = func
    return func


# --- Images -----------------------------------------------------------------

def to_pil(qimage: QImage) -> Image.Image:
    buffer = QBuffer()
    buffer.open(QIODevice.OpenModeFlag.WriteOnly)
    qimage.save(buffer, "PNG")
    return Image.open(io.BytesIO(bytes(buffer.data()))).convert("RGBA")


def grab(widget) -> Image.Image:
    """The widget as drawn, in physical pixels."""
    return to_pil(widget.grab().toImage())


def px(value):
    """Logical to physical pixels."""
    return int(round(value * SCALE))


def rounded_mask(size, radius, top_only=False):
    mask = Image.new("L", size, 0)
    ImageDraw.Draw(mask).rounded_rectangle(
        (0, 0, size[0] - 1, size[1] - 1), radius=radius, fill=255,
        corners=(True, True, not top_only, not top_only))
    return mask


def with_shadow(image: Image.Image, margin=SHADOW, offset=8, strength=90):
    """The image on a transparent canvas, with a soft shadow below it."""
    m, off = px(margin), px(offset)
    canvas = Image.new("RGBA", (image.width + 2 * m, image.height + 2 * m), (0, 0, 0, 0))
    alpha = Image.new("L", canvas.size, 0)
    alpha.paste(image.getchannel("A"), (m, m + off))
    alpha = alpha.filter(ImageFilter.GaussianBlur(px(margin) / 2.5))
    alpha = alpha.point(lambda a: a * strength // 255)
    shadow = Image.new("RGBA", canvas.size, (0, 0, 0, 255))
    shadow.putalpha(alpha)
    canvas.alpha_composite(shadow)
    canvas.alpha_composite(image, (m, m))
    return canvas


def title_bar(width, title, active=True):
    """A neutral title bar (not the one of a real system), in physical pixels."""
    image = QImage(width, px(TITLE_BAR), QImage.Format.Format_ARGB32)
    image.fill(QColor("#ececec" if active else "#f4f4f4"))
    painter = QPainter(image)
    font = QFont(QApplication.font())
    font.setPixelSize(px(13))
    font.setBold(True)
    painter.setFont(font)
    painter.setPen(QColor("#3a3a3a" if active else "#a0a0a0"))
    painter.drawText(QRect(0, 0, width, px(TITLE_BAR)), Qt.AlignmentFlag.AlignCenter, title)
    painter.setPen(QColor("#d0d0d0"))
    painter.drawLine(0, px(TITLE_BAR) - 1, width, px(TITLE_BAR) - 1)
    painter.end()
    return to_pil(image)


def framed(content: Image.Image, title: str, active=True, shadow=True):
    """A window: title bar, contents, rounded corners and shadow."""
    bar = title_bar(content.width, title, active)
    window = Image.new("RGBA", (content.width, bar.height + content.height))
    window.paste(bar, (0, 0))
    window.paste(content, (0, bar.height))
    window.putalpha(ImageChops.multiply(window.getchannel("A"),
                                        rounded_mask(window.size, px(RADIUS))))
    return with_shadow(window) if shadow else window


def popup(image: Image.Image):
    """A menu or a tooltip: rounded corners and a light shadow."""
    image = image.copy()
    image.putalpha(ImageChops.multiply(image.getchannel("A"),
                                       rounded_mask(image.size, px(5))))
    return with_shadow(image, margin=8, offset=2, strength=70)


def paste(base: Image.Image, overlay: Image.Image, at, shadow_margin=0):
    """Paste an overlay at a logical position of the base (shadow included)."""
    x, y = px(at[0]) - px(shadow_margin), px(at[1]) - px(shadow_margin)
    base.alpha_composite(overlay, (max(x, 0), max(y, 0)),
                         (max(-x, 0), max(-y, 0)))
    return base


def pointer(base: Image.Image, at):
    """Draw the mouse pointer with its tip at a logical position."""
    shape = [(0, 0), (0, 16), (4, 12.5), (6.8, 18.6), (9, 17.7), (6.3, 11.6), (11.5, 11.6)]
    x, y = px(at[0]), px(at[1])
    points = [(x + px(a), y + px(b)) for a, b in shape]
    draw = ImageDraw.Draw(base)
    draw.polygon(points, fill=(0, 0, 0, 255), outline=(255, 255, 255, 255), width=px(1))
    return base


def label(image: Image.Image, text, at, color="#e8a33a", size=30):
    """Draw an annotation centred at a logical position."""
    qimage = QImage(px(600), px(size * 1.6), QImage.Format.Format_ARGB32)
    qimage.fill(Qt.GlobalColor.transparent)
    painter = QPainter(qimage)
    font = QFont(QApplication.font())
    font.setPixelSize(px(size))
    font.setBold(True)
    painter.setFont(font)
    painter.setPen(QColor(color))
    painter.drawText(qimage.rect(), Qt.AlignmentFlag.AlignCenter, text)
    painter.end()
    overlay = to_pil(qimage)
    image.alpha_composite(overlay, (px(at[0]) - overlay.width // 2, px(at[1]) - overlay.height // 2))
    return image


def crop(image: Image.Image, box):
    """Crop a logical box (left, top, right, bottom)."""
    return image.crop(tuple(px(v) for v in box))


def data_file(name):
    return name.replace(" ", "_") + ".acs2"


def download(data_dir: Path):
    """Fetch the spectra that are not in the data folder yet."""
    data_dir.mkdir(parents=True, exist_ok=True)
    for name in SESSIONS:
        path = data_dir / data_file(name)
        if not path.exists():
            url = f"{DATA_RELEASE}/{data_file(name)}"
            print(f"Downloading {url}")
            urllib.request.urlretrieve(url, path)


# --- The GUI ----------------------------------------------------------------

class Studio:
    """Drives the GUI for the scenes."""

    def __init__(self, data_dir: Path):
        self.data_dir = data_dir
        self.app = QApplication.instance() or QApplication(sys.argv)
        # The default style of the offscreen platform draws no frame around
        # unchecked boxes; Fusion does, and is the same on every machine
        self.app.setStyle("Fusion")
        self.window = None
        self._analysed = None

    # Time and events
    def pump(self, seconds=0.3):
        end = time.time() + seconds
        while time.time() < end:
            self.app.processEvents()

    def wait_for_recipes(self, timeout=1800):
        end = time.time() + timeout
        while self.window.thread_pool.activeThreadCount() > 0:
            if time.time() > end:
                raise TimeoutError("recipe still running")
            self.pump(0.1)
        self.pump(0.5)

    # Windows and sessions
    def new_window(self, *sessions, size=WINDOW_SIZE):
        from astrocook.gui.main_window import MainWindowV2
        if self.window is not None:
            self.window.close()
            self.window.deleteLater()
            self.pump(0.2)
        # Each scene starts from the default state of the panels
        from PySide6.QtCore import QSettings
        QSettings("Astrocook", "AstrocookV2").remove("plot_controls")
        self.window = MainWindowV2(None, None)
        self.window.resize(*size)
        self.window.show()
        for name in sessions:
            self.open(name)
        return self.window

    def open(self, name):
        self.window.open_session_from_path(str(self.staged(name)))
        self.pump(1.0)

    def staged(self, name):
        """The file of a session, under the name the session must show."""
        if not hasattr(self, "_stage"):
            self._stage = Path(tempfile.mkdtemp(prefix="astrocook-shots-stage-"))
        target = self._stage / f"{name}.acs2"
        if not target.exists():
            source = self.data_dir / data_file(name)
            if not source.exists():
                raise FileNotFoundError(f"{source}: run with --download, "
                                        "see docs/screenshots/README.md")
            shutil.copyfile(source, target)
        return target

    def run(self, category, recipe, /, **params):
        """Run a recipe as its dialog does (parameters as strings)."""
        params = {k: str(v) for k, v in params.items()}
        self.window._on_recipe_requested(category, recipe, params, {})
        self.wait_for_recipes()

    def select_sessions(self, *names, current=None):
        """Select sessions in the list; the last (or ``current``) is shown."""
        from PySide6.QtCore import QItemSelectionModel
        view = self.window.session_list_view
        model = view.model()
        rows = [model.stringList().index(n) for n in names]
        current = rows[-1] if current is None else model.stringList().index(current)
        selection = view.selectionModel()
        view.clicked.emit(model.index(current, 0))
        self.pump(0.5)
        selection.clearSelection()
        for row in rows:
            selection.select(model.index(row, 0), QItemSelectionModel.SelectionFlag.Select)
        selection.setCurrentIndex(model.index(current, 0), QItemSelectionModel.SelectionFlag.NoUpdate)
        self.replot()

    def replot(self):
        self.window._trigger_replot()
        self.pump(1.0)
        # Cursor lines and continuum knots are blitted on mouse motion
        self.window.plot_viewer.canvas.draw_animated()
        self.pump(0.3)

    def limits(self, x=None, y=None):
        """Axis limits, as typed in the Axis Limits fields."""
        w = self.window
        if x is not None:
            w.xmin_input.setText(str(x[0])); w.xmax_input.setText(str(x[1]))
        if y is not None:
            w.ymin_input.setText(str(y[0])); w.ymax_input.setText(str(y[1]))
        w._on_set_custom_limits()
        self.pump(1.0)

    def open_panel(self, side):
        w = self.window
        card = w.right_sidebar_widget if side == "right" else w.left_sidebar_widget
        if card.width() == 0:
            (w._toggle_right_sidebar if side == "right" else w._toggle_left_sidebar)()
        self.pump(1.0)
        self.replot()

    def sections(self, **expanded):
        """Open or close sections of the plot controls (by title, _ for spaces)."""
        from astrocook.gui.main_window import CollapsibleSection
        found = {s.header.text(): s for s in
                 self.window.right_sidebar_widget.findChildren(CollapsibleSection)}
        for key, state in expanded.items():
            found[key.replace("_", " ")].setExpanded(state)
        self.pump(0.5)

    def recipe_dialog(self, category, name, values=None, focus=None):
        from astrocook.gui.recipe_dialog import RecipeDialog
        dialog = RecipeDialog(category, name, self.window.active_history.current_state, self.window)
        for key, value in (values or {}).items():
            widget = dialog.input_widgets[key]
            if hasattr(widget, "setText"):
                widget.setText(str(value))
            else:
                widget.setChecked(bool(value))
        dialog.show()
        if focus:
            dialog.input_widgets[focus].setFocus()
        self.pump(0.5)
        return dialog

    def context_menu(self, open_menu, active):
        """Open a context menu, grab it with an entry highlighted, close it.

        The GUI opens its menus with exec(), which waits for a choice: a
        timer grabs the menu from inside that wait and closes it. Returns
        the image and the geometry of the highlighted entry in the menu.
        """
        result = {}

        def grab_and_close():
            menu = QApplication.activePopupWidget()
            if not isinstance(menu, QMenu):
                QTimer.singleShot(100, grab_and_close)
                return
            action = next(a for a in menu.actions() if a.text().startswith(active))
            menu.setActiveAction(action)
            for _ in range(5):
                self.app.processEvents()
            result["image"] = grab(menu)
            result["item"] = menu.actionGeometry(action)
            menu.close()

        QTimer.singleShot(300, grab_and_close)
        open_menu()
        self.pump(0.3)
        return result["image"], result["item"]

    def analysed(self):
        """Open "My QSO" with continuum and fitted absorbers.

        The chain (properties, continuum, doublets, Lyman-alpha forest)
        takes about a minute and a half: it runs once, and later scenes open
        the saved session.
        """
        if self._analysed is None:
            self.new_window("My QSO")
            self.run("edit", "set_properties", name="My QSO", object="_current_",
                     z_em=Z_EM, resol=RESOL)
            self.run("continuum", "estimate_auto")
            self.run("absorbers", "doublets_auto")
            self.run("absorbers", "lya_auto")
            folder = Path(tempfile.mkdtemp(prefix="astrocook-shots-data-"))
            self._analysed = folder / "My QSO.acs2"
            self.window.active_history.current_state.save(str(self._analysed))
        self.new_window()
        self.window.open_session_from_path(str(self._analysed))
        self.pump(1.5)
        return self.window

    def at(self, widget, point=None):
        """Logical position of a point of a widget, in the main window."""
        point = widget.rect().center() if point is None else point
        p = widget.mapTo(widget.window(), point)
        return p.x(), p.y()

    def combo_popup(self, combo, active):
        """The list of a combo box, open, with an entry highlighted.

        Returns the image and its logical position in the window: the
        highlighted entry over the combo box, as Qt places it.
        """
        combo.showPopup()
        self.pump(0.3)
        view = combo.view()
        index = view.model().index(combo.findText(active), 0)
        view.setCurrentIndex(index)
        self.pump(0.2)
        container = view.window()
        image = grab(container)
        item = view.mapTo(container, view.visualRect(index).topLeft())
        combo.hidePopup()
        x, y = self.at(combo, QPoint(0, 0))
        return image, (x, y - item.y())

    def inspector(self, size=(1000, 620)):
        self.window._on_view_system_inspector()
        inspector = self.window.system_inspector
        inspector.resize(*size)
        self.pump(2.0)
        return inspector

    def select_components(self, inspector, *components):
        """Select the rows of the inspector of the given components."""
        from PySide6.QtCore import QItemSelectionModel
        model, view = inspector.proxy_model, inspector.table_view
        source = inspector.table_model

        def row_of(comp):
            for r in range(model.rowCount()):
                src = model.mapToSource(model.index(r, 0)).row()
                if source.get_component_at(src).uuid == comp.uuid:
                    return r
            raise LookupError(comp)
        rows = [row_of(c) for c in components]
        view.selectRow(rows[0])
        selection = view.selectionModel()
        for row in rows[1:]:
            selection.select(model.index(row, 0), QItemSelectionModel.SelectionFlag.Select
                             | QItemSelectionModel.SelectionFlag.Rows)
        view.scrollTo(model.index(rows[0], 0), view.ScrollHint.PositionAtCenter)
        self.pump(2.0)
        return [model.index(r, 0) for r in rows]

    def tooltip(self, text, widget):
        QToolTip.showText(widget.mapToGlobal(QPoint(0, 0)), text, widget)
        self.pump(0.3)
        label = next(w for w in self.app.topLevelWidgets()
                     if w.inherits("QTipLabel") and w.isVisible())
        image = grab(label)
        QToolTip.hideText()
        return image


# --- Scenes -----------------------------------------------------------------

@scene
def getting_started_main_view(s: Studio):
    w = s.new_window("My QSO")
    s.open_panel("left")
    s.open_panel("right")
    return framed(grab(w), w.windowTitle())


@scene
def getting_started_plot_controls(s: Studio):
    w = s.new_window("My QSO", size=(1000, 760))
    s.run("edit", "set_properties", name="My QSO", object="_current_",
          z_em=Z_EM, resol="_current_")
    s.limits(x=(409.5, 476.6))
    s.open_panel("right")
    s.sections(Redshift_Cursor=True)
    w.cursor_series_input.setText("Ly_a")
    w.cursor_z_input.setText("2.7310624")
    w.cursor_show_checkbox.setChecked(True)
    s.replot()
    image = framed(grab(w), w.windowTitle())
    return crop(image, (SHADOW + 1000 - 560, SHADOW + TITLE_BAR, 1000 + 2 * SHADOW,
                        SHADOW + TITLE_BAR + 760 + SHADOW))


@scene
def getting_started_session_list(s: Studio):
    w = s.new_window(*SESSIONS)
    s.open_panel("left")
    s.select_sessions("My other QSO")
    view = w.session_list_view
    row = view.visualRect(view.model().index(1, 0))
    menu, item = s.context_menu(lambda: w._on_session_list_context_menu(row.center()), "Close")
    x, y = s.at(view, row.center())
    image = grab(w)
    paste(image, popup(menu), (x + 60, y), shadow_margin=8)
    pointer(image, (x + 60 + item.right() - 30, y + item.center().y()))
    return crop(framed(image, w.windowTitle()), (0, 0, 560, 400))


@scene
def editing_multiselect_coadd(s: Studio):
    w = s.new_window(*SESSIONS)
    s.open_panel("left")
    s.select_sessions(*SESSIONS)
    s.limits(y=(-2000, 60000))
    view = w.session_list_view
    row = view.visualRect(view.model().index(1, 0))
    menu, item = s.context_menu(lambda: w._on_session_list_context_menu(row.center()), "Co-add")
    x, y = s.at(view, row.center())
    image = grab(w)
    paste(image, popup(menu), (x + 60, y), shadow_margin=8)
    pointer(image, (x + 60 + item.right() - 30, y + item.center().y()))
    return crop(framed(image, w.windowTitle()), (0, 0, 640, 520))


@scene
def getting_started_ghost_spectra(s: Studio):
    w = s.new_window(*SESSIONS)
    s.open_panel("left")
    s.select_sessions(*SESSIONS, current="My other QSO")
    s.limits(x=(393, 545), y=(-2000, 60000))
    return framed(grab(w), w.windowTitle())


@scene
def flux_smooth(s: Studio):
    w = s.new_window("My other QSO", size=(1000, 560))
    s.limits(x=(431.7, 442.4), y=(-2000, 20000))
    canvas = w.plot_viewer.canvas
    before = grab(canvas)
    s.run("flux", "smooth", sigma_kms=10.0)
    s.limits(x=(431.7, 442.4), y=(-2000, 20000))
    after = grab(canvas)
    # Left of 437 nm before the recipe, right of it after
    x_split = canvas.axes.transData.transform((437.0, 0))[0] / canvas.figure.bbox.width
    split = int(x_split * after.width)
    image = after.copy()
    image.paste(before.crop((0, 0, split, before.height)), (0, 0))
    dialog = s.recipe_dialog("flux", "smooth", {"sigma_kms": "10.0"}, focus="sigma_kms")
    field = dialog.input_widgets["sigma_kms"]
    fx, fy = s.at(field, field.rect().bottomLeft())
    tip = s.tooltip(field.toolTip() or "Standard deviation for Gaussian kernel (km/s)", field)
    framed_dialog = framed(grab(dialog), dialog.windowTitle())
    paste(image, framed_dialog, (190, 40), shadow_margin=SHADOW)
    tx, ty = 190 + fx + 60, 40 + TITLE_BAR + fy + 6
    paste(image, popup(tip), (tx, ty), shadow_margin=8)
    pointer(image, (tx - 4, ty - 14))
    dialog.close()
    label(image, "BEFORE  \u2192  AFTER", (canvas.width() * x_split, canvas.height() - 70))
    return image


@scene
def flux_running_std(s: Studio):
    w = s.new_window("My other QSO")
    s.run("flux", "calculate_running_std", input_col="y", output_col="running_std",
          window_pix=21)
    s.open_panel("right")
    w.aux_col_combo.setCurrentText("running_std")
    s.limits(x=(433.8, 439.6))
    s.replot()
    image = grab(w)
    popup_image, at = s.combo_popup(w.aux_col_combo, "running_std")
    paste(image, popup_image, at)
    combo_x, combo_y = s.at(w.aux_col_combo)
    pointer(image, (combo_x + 40, combo_y + 4))
    image = framed(image, w.windowTitle())
    return crop(image, (1000 + 2 * SHADOW - 560, 0, 1000 + 2 * SHADOW, 560))


@scene
def editing_set_properties(s: Studio):
    w = s.new_window("My QSO")
    s.limits(x=(440, 560))
    dialog = s.recipe_dialog("edit", "set_properties",
                             {"name": "My QSO", "object": "My QSO", "z_em": "2.730",
                              "resol": str(RESOL)},
                             focus="resol")
    field = dialog.input_widgets["resol"]
    tip = s.tooltip(field.toolTip(), field)
    image = grab(w.plot_viewer.canvas)
    fx, fy = s.at(field, field.rect().bottomLeft())
    paste(image, framed(grab(dialog), dialog.windowTitle()), (60, 90), shadow_margin=SHADOW)
    tx, ty = 60 + fx + 100, 90 + TITLE_BAR + fy + 6
    paste(image, popup(tip), (tx, ty), shadow_margin=8)
    pointer(image, (tx - 4, ty - 14))
    dialog.close()
    return image


@scene
def editing_apply_expression(s: Studio):
    w = s.new_window("My other QSO", "My QSO", size=(1000, 720))
    s.limits(x=(438, 462))
    dialog = s.recipe_dialog("edit", "apply_expression",
                             {"target_col": "F", "expression": "F / dF"}, focus="target_col")
    image = grab(w.plot_viewer.canvas)
    paste(image, framed(grab(dialog), dialog.windowTitle()), (120, 60), shadow_margin=SHADOW)
    dialog.close()
    return image


@scene
def continuum_auto(s: Studio):
    w = s.new_window("My QSO")
    s.run("edit", "set_properties", name="My QSO", object="_current_", z_em=Z_EM, resol=RESOL)
    s.limits(x=(429.5, 466))
    dialog = s.recipe_dialog("continuum", "estimate_auto", focus="smooth_len_lya")
    image = grab(w)
    at = (240, 160)
    paste(image, framed(grab(dialog), dialog.windowTitle()), at, shadow_margin=SHADOW)
    run = next(b for b in dialog.findChildren(QPushButton) if b.text() == "Run")
    rx, ry = s.at(run)
    pointer(image, (at[0] + rx + 20, at[1] + TITLE_BAR + ry))
    dialog.close()
    return framed(image, w.windowTitle(), active=False)


@scene
def absorbers_auto_pipeline(s: Studio):
    w = s.new_window("My QSO")
    s.run("edit", "set_properties", name="My QSO", object="_current_", z_em=Z_EM, resol=RESOL)
    s.run("continuum", "estimate_auto")
    s.limits(x=(429.5, 466))
    dialog = s.recipe_dialog("absorbers", "lya_auto", focus="score_threshold")
    image = grab(w)
    paste(image, framed(grab(dialog), dialog.windowTitle()), (180, 100), shadow_margin=SHADOW)
    pointer(image, (460, 470))
    dialog.close()
    return framed(image, w.windowTitle(), active=False)


@scene
def continuum_manual_edit(s: Studio):
    w = s.new_window("My QSO", size=(1000, 700))
    s.run("edit", "set_properties", name="My QSO", object="_current_", z_em=Z_EM, resol=RESOL)
    s.run("continuum", "estimate_auto")
    s.limits(x=(409.3, 472.9))
    s.open_panel("right")
    s.sections(Plot_Elements=False, Edit_Continuum=True)
    w._toggle_continuum_editor()
    s.pump(1.0)
    w.stride_entry.setText("2000")
    w._on_stride_text_changed()
    s.pump(1.0)
    s.replot()
    image = framed(grab(w), w.windowTitle())
    return crop(image, (1000 + 2 * SHADOW - 620, 200, 1000 + 2 * SHADOW, 700 + TITLE_BAR + 2 * SHADOW))


def _metal_system(session):
    """The strongest SiIV of the system at z ~ 2.3555 and its closest CIV."""
    comps = session.systs.components
    near = [c for c in comps if abs(c.z - 2.3555) < 0.002]
    siiv = max((c for c in near if c.series.startswith("SiIV")), key=lambda c: c.logN)
    civ = min((c for c in near if c.series.startswith("CIV")), key=lambda c: abs(c.z - siiv.z))
    return siiv, civ


def _inspector_scene(s: Studio, pair=False, menu_on=None, active=None):
    w = s.analysed()
    s.limits(x=(465, 472))
    siiv, civ = _metal_system(w.active_history.current_state)
    inspector = s.inspector()
    rows = s.select_components(inspector, *([siiv, civ] if pair else [siiv]))
    image = framed(grab(w), w.windowTitle(), active=False)
    insp = grab(inspector)
    if menu_on is not None:
        view = inspector.table_view
        col = [i for i in range(view.model().columnCount())
               if view.model().headerData(i, Qt.Orientation.Horizontal) == menu_on][0]
        rect = view.visualRect(rows[-1].siblingAtColumn(col))
        menu, item = s.context_menu(lambda: inspector._on_context_menu(rect.center()), active)
        x, y = s.at(view, rect.center())
        paste(insp, popup(menu), (x, y), shadow_margin=8)
        pointer(insp, (x + item.right() - 30, y + item.center().y()))
    else:
        _component_tooltip(s, inspector, insp)
    inspector.hide()
    inspector_framed = framed(insp, inspector.windowTitle())
    canvas = Image.new("RGBA", (image.width + px(140), image.height + px(30)), (0, 0, 0, 0))
    canvas.alpha_composite(image, (0, 0))
    canvas.alpha_composite(inspector_framed, (px(140), px(30)))
    return canvas


def _component_tooltip(s: Studio, inspector, image):
    """The tooltip of a tick of the velocity plot, as on mouse over."""
    from matplotlib.backend_bases import MouseEvent
    plot = inspector.vel_plot
    shown = []
    original = QToolTip.showText
    QToolTip.showText = staticmethod(lambda pos, text, *a: shown.append(text))
    try:
        ax, v, comp = next(t for t in plot._tick_map_visuals if t[1] > 100) \
            if any(t[1] > 100 for t in plot._tick_map_visuals) else plot._tick_map_visuals[-1]
        x, y = ax.transData.transform((v, 0.5))
        event = MouseEvent("motion_notify_event", plot.canvas, x, y)
        plot.on_mouse_move(event)
    finally:
        QToolTip.showText = original
    if not shown:
        return
    canvas = plot.canvas
    ratio = canvas.devicePixelRatioF()
    cx = x / ratio
    cy = canvas.height() - y / ratio
    tip = s.tooltip(shown[-1], canvas)
    px_, py_ = s.at(canvas, QPoint(int(cx), int(cy)))
    paste(image, popup(tip), (px_ + 8, py_ + 12), shadow_margin=8)
    pointer(image, (px_, py_))


@scene
def absorbers_system_inspector(s: Studio):
    return _inspector_scene(s)


@scene
def absorbers_linking(s: Studio):
    return _inspector_scene(s, pair=True, menu_on="b", active="Link b to SiIV (Thermal)")


# --- Main -------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA,
                        help="folder with 'My QSO.acs2' and 'My other QSO.acs2'")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output folder")
    parser.add_argument("--only", nargs="+", metavar="SCENE", help="scenes to make")
    parser.add_argument("--list", action="store_true", help="list the scenes")
    parser.add_argument("--download", action="store_true",
                        help="download the spectra missing from --data")
    args = parser.parse_args(argv)
    if args.list:
        print("\n".join(SCENES))
        return 0
    names = args.only or list(SCENES)
    unknown = [n for n in names if n not in SCENES]
    if unknown:
        parser.error(f"unknown scenes: {', '.join(unknown)}")
    if args.download:
        download(args.data)
    studio = Studio(args.data)
    args.out.mkdir(parents=True, exist_ok=True)
    for name in names:
        start = time.time()
        image = SCENES[name](studio)
        path = args.out / f"{name}.png"
        image.save(path, optimize=True)
        print(f"{path} ({image.width}x{image.height}, {time.time() - start:.0f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
