"""Unit tests for command line script functions."""

import importlib.resources

import numpy as np
import pytest

from snf_simulations.cask import Cask
from snf_simulations.physics import calculate_event_rate, calculate_flux_at_distance
from snf_simulations.spec import Spectrum

from . import data

# Suppress warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers


def _load_proportions(reactor: str) -> dict[str, float]:
    """Load in example reactor isotope proportions."""
    with importlib.resources.path(
        data, f"{reactor.capitalize()}_proportions.csv"
    ) as filepath:
        reactor_data = np.genfromtxt(
            filepath,
            delimiter=",",
            skip_header=1,
            dtype=str,
        )
    return {str(d[0]): float(d[1]) for d in reactor_data}


def _load_output(filename: str) -> tuple[np.ndarray, np.ndarray]:
    """Load data from a CSV output file."""
    with importlib.resources.path(data, filename) as path:
        spec_data = np.loadtxt(
            path,
            delimiter=",",
            skiprows=1,
        )
    energy = spec_data[:, 0].tolist()
    flux = spec_data[:, 1].tolist()
    return energy, flux


@pytest.mark.parametrize("reactor", ["sizewell", "hartlepool"])
def test_single_cask(reactor: str) -> None:
    """Test single cask spectrum output."""
    cask_mass = 10000  # 10 tonnes

    # Create the Cask for the given reactor.
    # Note the proportions in the file assumed a 1 tonne cask mass,
    # not the more accurate 1.135 tonne total.
    # Also they're from the 24 hour simulation time, not the lowest.
    # However we didn't take that into account when creating the reference data,
    # so we still set the initial_cooling_time to zero.
    isotope_proportions = _load_proportions(reactor)
    isotope_masses = {
        isotope: proportion * cask_mass
        for isotope, proportion in isotope_proportions.items()
    }
    initial_cooling_time = 0  # 24 * _UNITS_TO_YEARS["HOURS"]
    cask = Cask(isotope_masses, initial_cooling_time, name=f"{reactor}_cask")

    # Get the Spectra for 6 months after removal from the core.
    spec = cask.get_total_spectrum(cooling_time=0.5)

    # Compare to the reference data file
    with importlib.resources.path(data, f"{reactor.capitalize()}_single.csv") as path:
        spec_ref = Spectrum.from_file(path)

    assert len(spec.energy) == len(spec_ref.energy), (
        f"Energy data has wrong length ({len(spec.energy)} vs {len(spec_ref.energy)})"
    )
    assert len(spec.flux) == len(spec_ref.flux), (
        f"Flux data has wrong length ({len(spec.flux)} vs {len(spec_ref.flux)})"
    )
    assert np.allclose(spec.energy, spec_ref.energy), "Reference energy does not match"
    assert np.allclose(spec.flux, spec_ref.flux), "Reference flux does not match"

    # Calculate the single cask flux at 40m
    # Returns a single value of flux cm^-2 s^-1
    total_flux = spec.integrate(1806, 6000)
    flux_at_40m = calculate_flux_at_distance(total_flux, distance=40)
    assert isinstance(flux_at_40m, float), "Single flux is not a float"
    if reactor == "sizewell":
        flux_ref = 1.191153e09
    elif reactor == "hartlepool":
        flux_ref = 4.909233e08
    assert np.isclose(flux_at_40m, flux_ref), (
        f"Single flux value does not match: ({flux_ref:e} vs {flux_at_40m:e})"
    )

    # Calculate event rates in a detector at 40m using the flux spectrum
    rate_lower, rate_upper = calculate_event_rate(flux_at_40m, 0.2, 0.4)
    assert isinstance(rate_lower, float), "Lower event rate is not a float"
    assert isinstance(rate_upper, float), "Upper event rate is not a float"
    if reactor == "sizewell":
        rate_lower_ref = 1.315033e-07
        rate_upper_ref = 2.630066e-07
    elif reactor == "hartlepool":
        rate_lower_ref = 5.419793e-08
        rate_upper_ref = 1.083959e-07
    assert np.isclose(rate_lower, rate_lower_ref, atol=1e-15), (
        f"Lower event rate does not match: ({rate_lower_ref:e} vs {rate_lower:e})"
    )
    assert np.isclose(rate_upper, rate_upper_ref, atol=1e-15), (
        f"Upper event rate does not match: ({rate_upper_ref:e} vs {rate_upper:e})"
    )


@pytest.mark.parametrize("reactor", ["sizewell", "hartlepool"])
def test_multiple_casks(reactor: str) -> None:
    """Test multiple cask spectrum output."""
    cask_mass = 10000
    n_casks = 10
    if reactor == "sizewell":
        cooling_times = [0.5, 5, 10, 20]
    elif reactor == "hartlepool":
        cooling_times = [3, 7, 15, 19]

    # Create a combined total Cask for the given reactor (see notes above).
    isotope_proportions = _load_proportions(reactor)
    isotope_masses = {
        isotope: proportion * cask_mass * n_casks
        for isotope, proportion in isotope_proportions.items()
    }
    initial_cooling_time = 0  # 24 * _UNITS_TO_YEARS["HOURS"]
    cask = Cask(isotope_masses, initial_cooling_time, name=f"{reactor}_cask")

    # Create the Spectra for each set of casks at the specified times.
    spectra = []
    for cooling_time in cooling_times:
        spec = cask.get_total_spectrum(cooling_time=cooling_time)
        spectra.append(spec)

    # Combine the spectra from all casks to get the total spectrum for all 40 casks.
    spec_multiple = spectra[0]
    for spec in spectra[1:]:
        spec_multiple = spec_multiple + spec

    # Compare to the reference data file
    with importlib.resources.path(data, f"{reactor.capitalize()}_multiple.csv") as path:
        spec_ref = Spectrum.from_file(path)
    assert len(spec_ref.energy) == len(spec_multiple.energy), (
        f"Energy data has wrong length "
        f"({len(spec_ref.energy)} vs {len(spec_multiple.energy)})"
    )
    assert len(spec_ref.flux) == len(spec_multiple.flux), (
        f"Flux data has wrong length "
        f"({len(spec_ref.flux)} vs {len(spec_multiple.flux)})"
    )
    assert np.allclose(spec_ref.energy, spec_multiple.energy), (
        "Reference energy does not match"
    )
    assert np.allclose(spec_ref.flux, spec_multiple.flux), (
        "Reference flux does not match"
    )

    # Calculate the total cask flux at 40m
    total_flux = spec_multiple.integrate(1806, 6000)
    flux_at_40m = calculate_flux_at_distance(total_flux, distance=40)
    assert isinstance(flux_at_40m, float), "Total flux is not a float"
    if reactor == "sizewell":
        flux_ref = 1.297917e10
    elif reactor == "hartlepool":
        flux_ref = 1.431710e09
    assert np.isclose(flux_at_40m, flux_ref), (
        f"Total flux value does not match: ({flux_ref:e} vs {flux_at_40m:e})"
    )

    # Calculate event rates in a detector at 40m using the flux spectrum
    rate_lower, rate_upper = calculate_event_rate(flux_at_40m, 0.2, 0.4)
    assert isinstance(rate_lower, float), "Lower event rate is not a float"
    assert isinstance(rate_upper, float), "Upper event rate is not a float"
    if reactor == "sizewell":
        rate_lower_ref = 1.432901e-06
        rate_upper_ref = 2.865801e-06
    elif reactor == "hartlepool":
        rate_lower_ref = 1.580608e-07
        rate_upper_ref = 3.161216e-07
    assert np.isclose(rate_lower, rate_lower_ref, atol=1e-15), (
        f"Lower event rate does not match: ({rate_lower_ref:e} vs {rate_lower:e})"
    )
    assert np.isclose(rate_upper, rate_upper_ref, atol=1e-15), (
        f"Upper event rate does not match: ({rate_upper_ref:e} vs {rate_upper:e})"
    )
