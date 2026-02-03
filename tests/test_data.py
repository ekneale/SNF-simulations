"""Unit tests for data loading functions."""

import numpy as np
import pytest

from snf_simulations.data import (
    get_reactors,
    load_antineutrino_data,
    load_isotope_data,
    load_reactor_data,
    load_spectrum,
)

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts


def test_get_reactors():
    """Test that reactor list can be retrieved."""
    reactors = get_reactors()
    assert isinstance(reactors, list), "Reactors should be a list"
    assert len(reactors) > 0, "Should have at least one reactor"
    assert all(isinstance(r, str) for r in reactors), "All reactors should be strings"
    assert "sizewell" in reactors, "Sizewell should be in reactor list"
    assert "hartlepool" in reactors, "Hartlepool should be in reactor list"


def test_load_reactor_data_sizewell():
    """Test loading reactor data for Sizewell."""
    data = load_reactor_data("sizewell")
    assert isinstance(data, dict), "Loaded data should be a dictionary"
    assert len(data) > 0, "Should have loaded isotope data"
    assert all(isinstance(k, str) for k in data), "All keys should be isotope names"
    assert all(isinstance(v, float) for v in data.values()), (
        "All values should be floats"
    )
    assert all(v >= 0 for v in data.values()), "All proportions should be non-negative"


def test_load_reactor_data_hartlepool():
    """Test loading reactor data for Hartlepool."""
    data = load_reactor_data("hartlepool")
    assert isinstance(data, dict), "Loaded data should be a dictionary"
    assert len(data) > 0, "Should have loaded isotope data"


def test_load_reactor_data_invalid():
    """Test that loading invalid reactor raises ValueError."""
    with pytest.raises(ValueError, match=r"Reactor.*data file not found"):
        load_reactor_data("invalid_reactor")


def test_load_isotope_data_all():
    """Test loading all isotope data."""
    molar_masses, half_lives = load_isotope_data()
    assert isinstance(molar_masses, dict), "Molar masses should be a dictionary"
    assert isinstance(half_lives, dict), "Half lives should be a dictionary"
    assert len(molar_masses) > 0, "Should have loaded molar masses"
    assert len(half_lives) > 0, "Should have loaded half lives"
    assert len(molar_masses) == len(half_lives), (
        "Molar masses and half lives should have same length"
    )

    # Check data types and values
    assert all(isinstance(k, str) for k in molar_masses), (
        "All molar mass keys should be isotope names"
    )
    assert all(isinstance(v, int) for v in molar_masses.values()), (
        "All molar masses should be integers"
    )
    assert all(v > 0 for v in molar_masses.values()), (
        "All molar masses should be positive"
    )

    assert all(isinstance(k, str) for k in half_lives), (
        "All half life keys should be isotope names"
    )
    assert all(isinstance(v, float) for v in half_lives.values()), (
        "All half lives should be floats"
    )
    assert all(v > 0 for v in half_lives.values()), "All half lives should be positive"


def test_load_isotope_data():
    """Test loading specific isotopes."""
    isotopes = ["Sr90", "Y90"]
    molar_masses, half_lives = load_isotope_data(isotopes)
    assert len(molar_masses) == 2, "Should have loaded 2 isotopes"
    assert len(half_lives) == 2, "Should have loaded 2 isotopes"
    assert "Sr90" in molar_masses, "Sr90 should be in molar masses"
    assert "Y90" in molar_masses, "Y90 should be in molar masses"
    assert "Sr90" in half_lives, "Sr90 should be in half lives"
    assert "Y90" in half_lives, "Y90 should be in half lives"
    assert molar_masses["Sr90"] == 90, "Molar mass of Sr90 should be 90 g/mol"
    assert molar_masses["Y90"] == 90, "Molar mass of Y90 should be 90 g/mol"
    assert half_lives["Sr90"] == 28.91, (
        "Half life of Sr90 should be approximately 28.91 years"
    )
    assert half_lives["Y90"] == 0.0073, (
        "Half life of Y90 should be approximately 0.0073 years"
    )


def test_load_spectrum():
    """Test loading spectrum data for a valid isotope."""
    spectrum = load_spectrum("Sr90")
    assert isinstance(spectrum, np.ndarray), "Spectrum should be a numpy array"
    assert spectrum.shape[1] == 3, (
        "Spectrum should have 3 columns (energy, dN/dE, uncertainty)"
    )
    assert spectrum.shape[0] > 0, "Spectrum should have at least one row"

    # Check that energy values are positive and increasing
    energy = spectrum[:, 0]
    assert np.all(energy >= 0), "Energy values should be positive"
    assert np.all(np.diff(energy) > 0), "Energy should be increasing"

    # Check that flux values are non-negative
    flux = spectrum[:, 1]
    assert np.all(flux >= 0), "Flux values should be non-negative"

    # Check that uncertainties are non-negative
    uncertainty = spectrum[:, 2]
    assert np.all(uncertainty >= 0), "Uncertainty values should be non-negative"


def test_load_spectrum_invalid():
    """Test that loading invalid isotope raises ValueError."""
    with pytest.raises(ValueError, match="Spectrum data file.*not found"):
        load_spectrum("InvalidIsotope")


def test_load_antineutrino_data():
    """Test loading antineutrino data for multiple isotopes."""
    isotopes = ["Sr90", "Y90"]
    data = load_antineutrino_data(isotopes)
    assert isinstance(data, dict), "Loaded data should be a dictionary"
    assert len(data) == 2, "Should have loaded 2 isotopes"

    for isotope in isotopes:
        assert isotope in data, f"{isotope} should be in loaded data"
        spectrum = data[isotope]
        assert isinstance(spectrum, np.ndarray), (
            f"{isotope} spectrum should be numpy array"
        )
        assert spectrum.shape[1] == 3, f"{isotope} spectrum should have 3 columns"
