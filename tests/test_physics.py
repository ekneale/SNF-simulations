"""Unit tests for physics calculations."""

import numpy as np
import ROOT

from snf_simulations.physics import (
    DecayChain,
    calculate_event_rate,
    calculate_flux,
    get_decay_mass,
)

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts


def test_decay_chain_defaults():
    """Test DecayChain default branching ratio."""
    chain = DecayChain("U235", "Th231")
    assert chain.parent == "U235", "Parent isotope should be U235"
    assert chain.daughter == "Th231", "Daughter isotope should be Th231"
    assert chain.branching_ratio == 1.0, "Default branching ratio should be 1.0"


def test_get_decay_mass_zero_time():
    """Test that decay mass is zero at time zero."""
    daughter_mass = get_decay_mass(
        time_elapsed=0,
        parent_mass=10.0,
        parent_half_life=2.0,
        daughter_half_life=3.0,
        branching_ratio=1.0,
    )
    assert np.isclose(daughter_mass, 0.0), "Daughter mass should be zero at time zero"


def test_get_decay_mass_positive_time():
    """Test decay mass against a manually computed value."""
    time_elapsed = 5.0
    parent_mass = 10.0
    parent_half_life = 2.0
    daughter_half_life = 3.0
    branching_ratio = 0.5
    daughter_mass = get_decay_mass(
        time_elapsed=time_elapsed,
        parent_mass=parent_mass,
        parent_half_life=parent_half_life,
        daughter_half_life=daughter_half_life,
        branching_ratio=branching_ratio,
    )

    # Bateman equation for daughter isotope mass
    parent_decay_constant = np.log(2) / parent_half_life
    daughter_decay_constant = np.log(2) / daughter_half_life
    daughter_mass_ref = (
        branching_ratio
        * (parent_decay_constant / (daughter_decay_constant - parent_decay_constant))
        * parent_mass
        * (
            np.exp(-parent_decay_constant * time_elapsed)
            - np.exp(-daughter_decay_constant * time_elapsed)
        )
    )

    assert np.isclose(daughter_mass, daughter_mass_ref), (
        "Daughter mass does not match reference value"
    )


def test_calculate_flux():
    """Test flux calculation against expected value."""
    integral_value = 8.0e9
    spec = ROOT.TH1D("spec", "", 6000, 0, 6000)
    spec.SetBinContent(2000, integral_value)
    distance = 40.0
    flux = calculate_flux(spec, distance)
    flux_ref = integral_value / (4 * np.pi * (distance * 100) ** 2)
    assert np.isclose(flux, flux_ref), "Calculated flux does not match reference value"


def test_calculate_event_rate():
    """Test event rate calculation with default efficiency bounds."""
    flux = 2.0e9
    rate_lower, rate_upper = calculate_event_rate(flux)

    # The detector properties are currently hardcoded in the function
    detector_volume = 1.52 * 1.52 * 0.7 * 1e6
    number_of_protons = detector_volume * 4.6e22
    event_rate = number_of_protons * 1e-44 * flux

    assert np.isclose(rate_lower, event_rate * 0.2), "Lower event rate does not match"
    assert np.isclose(rate_upper, event_rate * 0.4), "Upper event rate does not match"
