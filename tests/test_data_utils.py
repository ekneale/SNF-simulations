"""Unit tests for data utility functions."""

import pytest

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
