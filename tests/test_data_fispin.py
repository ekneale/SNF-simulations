"""Unit tests for FISPIN data functions."""

from pathlib import Path

import pytest

from snf_simulations.data.fispin import (
    _convert_sim_time_to_years,
    get_isotope_masses,
    load_tabqfile,
)

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers


EXAMPLE_TABQ_CONTENT = (
    "*** TIME    1.200E+01 HOURS\n"
    "ALL-NUC      GRAMS\n"
    "SR90         3.0\n"
    "CS137        3.0\n"
    "TOTAL        6.0\n"
    "*** TIME    2.000E+00 DAYS\n"
    "ALL-NUC      GRAMS\n"
    "SR90         2.0\n"
    "CS137        8.0\n"
    "TOTAL        10.0\n"
)


def _write_tabqfile(tmp_path: Path) -> Path:
    """Create a minimal synthetic .tbQ file with two time sections."""
    filepath = tmp_path / "sample.tbQ"
    filepath.write_text(EXAMPLE_TABQ_CONTENT, encoding="utf-8")
    return filepath


def test_convert_sim_time_to_years() -> None:
    """Test conversion of supported simulation time units to years."""
    assert _convert_sim_time_to_years("1.000E+00 YEARS") == pytest.approx(1.0)
    assert _convert_sim_time_to_years("3.652425E+02 DAYS") == pytest.approx(1.0)
    assert _convert_sim_time_to_years("8.76582E+03 HOURS") == pytest.approx(1.0)


def test_load_tabqfile(tmp_path: Path) -> None:
    """Test .tbQ loading returns dataframes for each time section."""
    filepath = _write_tabqfile(tmp_path)

    time_dfs = load_tabqfile(filepath)

    assert set(time_dfs) == {"2.000E+00 DAYS", "1.200E+01 HOURS"}

    df = time_dfs["2.000E+00 DAYS"]
    assert list(df.columns) == ["ALL-NUC", "GRAMS"]
    assert df["ALL-NUC"].tolist() == ["Sr90", "Cs137"]
    assert df["GRAMS"].tolist() == pytest.approx([2.0, 8.0])


def test_load_tabqfile_from_string() -> None:
    """Test .tbQ loading from a string input."""
    time_dfs = load_tabqfile(EXAMPLE_TABQ_CONTENT)

    assert set(time_dfs) == {"2.000E+00 DAYS", "1.200E+01 HOURS"}

    df = time_dfs["1.200E+01 HOURS"]
    assert list(df.columns) == ["ALL-NUC", "GRAMS"]
    assert df["ALL-NUC"].tolist() == ["Sr90", "Cs137"]
    assert df["GRAMS"].tolist() == pytest.approx([3.0, 3.0])


def test_load_tabqfile_invalid_path() -> None:
    """Test that loading from an invalid file path raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_tabqfile(Path("non_existent_file.tbQ"))


def test_get_isotope_masses_default_time(tmp_path: Path) -> None:
    """Test default behavior selects earliest simulation time in the file."""
    filepath = _write_tabqfile(tmp_path)

    masses, cooling_time = get_isotope_masses(filepath)

    assert cooling_time == pytest.approx(0.5 / 365.2425)
    assert masses["Sr90"] == pytest.approx(0.003)
    assert masses["Cs137"] == pytest.approx(0.003)
    assert sum(masses.values()) == pytest.approx(0.006)


def test_get_isotope_masses_selected_time(tmp_path: Path) -> None:
    """Test requesting a specific time returns matching isotope masses."""
    filepath = _write_tabqfile(tmp_path)

    masses, cooling_time = get_isotope_masses(filepath, "2.000E+00 DAYS")

    assert cooling_time == pytest.approx(2.0 / 365.2425)
    assert masses["Sr90"] == pytest.approx(0.002)
    assert masses["Cs137"] == pytest.approx(0.008)


def test_get_isotope_masses_invalid_time(tmp_path: Path) -> None:
    """Test requesting a missing simulation time raises ValueError."""
    filepath = _write_tabqfile(tmp_path)

    with pytest.raises(ValueError, match=r"Specified time string"):
        get_isotope_masses(filepath, "1.000E+01 YEARS")
