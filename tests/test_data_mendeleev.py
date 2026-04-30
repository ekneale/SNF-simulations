"""Unit tests for isotope data using the mendeleev package."""

import pytest

from snf_simulations.data.mendeleev import get_isotope_properties

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers


def test_get_isotope_properties(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test loading isotope data."""

    class _MockIsotope:
        mass = 90.0
        half_life = 100.0
        half_life_unit = "year"

    monkeypatch.setattr(
        "snf_simulations.data.mendeleev.isotope", lambda *_: _MockIsotope()
    )

    isotope_properties = get_isotope_properties("Y90")
    assert isotope_properties["molar_mass"] == 90.0, "Molar mass should be 90 g/mol"
    assert isotope_properties["half_life"] == 100.0, "Half life should be 100 years"


def test_get_isotope_properties_real() -> None:
    """Test loading specific isotope data."""
    # We can't test every possible isotope, but we can do one to check there
    # aren't any obvious problems.
    isotope_properties = get_isotope_properties("Y90")
    assert isotope_properties["molar_mass"] == pytest.approx(89.9, rel=1e-3), (
        "Molar mass of Y90 should be approximately 89.9 g/mol"
    )
    assert isotope_properties["half_life"] == pytest.approx(0.0073, rel=1e-3), (
        "Half life of Y90 should be approximately 0.0073 years"
    )


def test_get_isotope_properties_converts_to_years(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test half-life conversion from non-year units into years."""

    class _MockIsotope:
        mass = 90.0
        half_life = 24.0
        half_life_unit = "hour"

    monkeypatch.setattr(
        "snf_simulations.data.mendeleev.isotope", lambda *_: _MockIsotope()
    )

    isotope_properties = get_isotope_properties("Y90")

    assert isotope_properties["molar_mass"] == 90.0, "Molar mass should be 90 g/mol"
    assert isotope_properties["half_life"] == pytest.approx(
        24.0 * 3600.0 / 31_556_926.0
    ), "Half life in years should be 24 hours converted to years"


def test_get_isotope_properties_unsupported_unit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test unsupported half-life unit raises ValueError."""

    class _MockIsotope:
        mass = 1.0
        half_life = 1.0
        half_life_unit = "fortnight"

    monkeypatch.setattr(
        "snf_simulations.data.mendeleev.isotope", lambda *_: _MockIsotope()
    )

    with pytest.raises(ValueError, match=r"Unsupported half-life unit"):
        get_isotope_properties("Y90")
