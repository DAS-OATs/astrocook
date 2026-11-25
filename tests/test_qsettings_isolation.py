"""The tests never touch the user's real settings (phase 7-quater, U2)."""
from astrocook.gui.settings import app_settings


def test_settings_go_to_a_folder_of_the_test(tmp_path):
    settings = app_settings()
    settings.setValue("probe", 1)
    settings.sync()
    assert str(tmp_path) in settings.fileName()


def test_settings_do_not_leak_between_tests_1():
    app_settings().setValue("leak", "yes")


def test_settings_do_not_leak_between_tests_2():
    assert app_settings().value("leak") is None
