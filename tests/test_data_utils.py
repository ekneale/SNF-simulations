"""Unit tests for data utility functions."""

from pathlib import Path

import pytest

from snf_simulations.data import get_example_tbq_path
from snf_simulations.data.utils import _parse_isotope

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers


def test_parse_isotope() -> None:
    """Test parsing isotope names in both supported orders."""
    assert _parse_isotope("Ru106") == ("Ru", 106)
    assert _parse_isotope("106Ru") == ("Ru", 106)
    assert _parse_isotope("y90") == ("y", 90)

    with pytest.raises(ValueError, match=r"Isotope format not recognized"):
        _parse_isotope("Ru-106")
    with pytest.raises(ValueError, match=r"Isotope format not recognized"):
        _parse_isotope("Ru")
    with pytest.raises(ValueError, match=r"Isotope format not recognized"):
        _parse_isotope("106")


def test_get_example_tbq_path() -> None:
    """Test that get_example_tbq_path returns expected Path."""
    path = get_example_tbq_path()
    assert isinstance(path, Path)
    assert path.exists()
    assert path.is_file()
    assert "snf_simulations" in path.parts
    assert "data" in path.parts
    assert path.name == "example.tbQ"
