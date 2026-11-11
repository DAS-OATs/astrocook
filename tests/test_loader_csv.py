"""The loader of CSV files (F97, F100, after phase 6).

The first three columns are x (nm), y and dy, taken by position, with or
without a header. Before: the first line was always the header, and the
columns were asked for by the names col1, col2, col3 (F100); a file of
fewer than 10 lines went to the loader of FITS files (F97).
"""
import numpy as np
import pytest

from astrocook.session.session import load_session_from_file
from astrocook.io.loaders import load_simple_csv

ROWS = "500.0,1.0,0.1\n501.0,1.5,0.2\n502.0,2.0,0.3\n"


def _values(data):
    return (list(data.x.values), list(data.y.values), list(data.dy.values))


EXPECTED = ([500.0, 501.0, 502.0], [1.0, 1.5, 2.0], [0.1, 0.2, 0.3])


@pytest.mark.parametrize("header", ["col1,col2,col3\n", "wave,flux,err\n", "",
                                    "# a comment\n", "# a comment\nwave,flux,err\n"])
def test_csv_with_or_without_a_header(tmp_path, header):
    path = tmp_path / "s.csv"
    path.write_text(header + ROWS)
    assert _values(load_simple_csv(str(path))) == EXPECTED


def test_more_columns_take_the_first_three(tmp_path):
    path = tmp_path / "s.csv"
    path.write_text("wave,flux,err,cont\n" + ROWS.replace("\n", ",9.0\n"))
    assert _values(load_simple_csv(str(path))) == EXPECTED


def test_fewer_than_three_columns(tmp_path):
    path = tmp_path / "s.csv"
    path.write_text("500.0,1.0\n501.0,1.5\n")
    with pytest.raises(ValueError, match="2 columns: three are needed"):
        load_simple_csv(str(path))


def test_csv_without_a_header_opens_as_a_session(tmp_path):
    # Short (F97) and without a header (F100)
    path = tmp_path / "s.csv"
    path.write_text(ROWS)
    session = load_session_from_file(str(path), "s", None, "auto")
    assert session != 0
    np.testing.assert_array_equal(session.spec.x.value, EXPECTED[0])
    assert str(session.spec.x.unit) == "nm"
