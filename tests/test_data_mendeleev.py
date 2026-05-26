"""Unit tests for isotope data using the mendeleev package."""

from collections.abc import Iterator

import numpy as np
import pytest

from snf_simulations.data.mendeleev import (
    _get_isotope_properties_cached,
    get_isotope_properties,
)

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers


@pytest.fixture(autouse=True)
def _clear_isotope_properties_cache() -> Iterator[None]:
    """Reset the mendeleev properties cache before and after each test."""
    _get_isotope_properties_cached.cache_clear()
    yield  # The test runs here
    _get_isotope_properties_cached.cache_clear()


def test_get_isotope_properties(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test loading isotope data."""

    class _MockIsotope:
        mass = 90.0
        half_life = 100.0
        half_life_unit = "year"
        decay_modes = [type("DecayMode", (), {"mode": "B-"})()]

    monkeypatch.setattr(
        "snf_simulations.data.mendeleev.isotope", lambda *_: _MockIsotope()
    )

    isotope_properties = get_isotope_properties("Y90")
    assert isotope_properties["molar_mass"] == 90.0, "Molar mass should be 90 g/mol"
    assert isotope_properties["half_life"] == 100.0, "Half life should be 100 years"
    assert isotope_properties["decay_modes"] == ["B-"], "Decay modes should be ['B-']"


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
    assert isotope_properties["decay_modes"] == ["B-"], (
        "Decay modes of Y90 should be ['B-']"
    )

    # Also test an isotope with a stable half-life to check that np.inf is returned
    isotope_properties = get_isotope_properties("H1")
    assert isotope_properties["molar_mass"] == pytest.approx(1.00784, rel=1e-3), (
        "Molar mass of H1 should be approximately 1.00784 g/mol"
    )
    assert isotope_properties["half_life"] == np.inf, (
        "Half life of stable isotope should be np.inf"
    )
    assert isotope_properties["decay_modes"] == [], (
        "Decay modes of stable isotope should be empty"
    )


def test_get_isotope_properties_converts_to_years(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test half-life conversion from non-year units into years."""

    class _MockIsotope:
        mass = 90.0
        half_life = 24.0
        half_life_unit = "hour"
        decay_modes = [type("DecayMode", (), {"mode": "B-"})()]

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
        decay_modes = [type("DecayMode", (), {"mode": "B-"})()]

    monkeypatch.setattr(
        "snf_simulations.data.mendeleev.isotope", lambda *_: _MockIsotope()
    )

    with pytest.raises(ValueError, match=r"Unsupported half-life unit"):
        get_isotope_properties("Y90")


def test_get_isotope_properties_cache(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test that repeated lookups should reuse results from the cache."""
    call_count = 0

    class _MockIsotope:
        mass = 90.0
        half_life = 100.0
        half_life_unit = "year"
        decay_modes = [type("DecayMode", (), {"mode": "B-"})()]

    def _mock_isotope(*_: object) -> _MockIsotope:
        # By using a nonlocal variable we can check how many times this
        # function is called, which should be only once due to caching.
        nonlocal call_count
        call_count += 1
        return _MockIsotope()

    monkeypatch.setattr("snf_simulations.data.mendeleev.isotope", _mock_isotope)

    # First try to get the properties for the same isotope twice
    first = get_isotope_properties("Y90")
    second = get_isotope_properties("Y90")
    assert first == second
    assert call_count == 1, (
        "mendeleev.isotope should be called once for a cached isotope"
    )

    # Now try after cleaning the cache in between
    _get_isotope_properties_cached.cache_clear()
    third = get_isotope_properties("Y90")
    assert third == first
    assert call_count == 2, (
        "mendeleev.isotope should be called again after cache is cleared"
    )

    # Finally, try again without clearing the cache to confirm it is still cached
    fourth = get_isotope_properties("Y90")
    assert fourth == first
    assert call_count == 2, (
        "mendeleev.isotope should not be called again for a cached isotope"
    )
