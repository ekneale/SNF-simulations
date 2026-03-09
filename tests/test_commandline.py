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
    cask_mass = 10000  # 10 tonne casks
    removal_times = [0, 0.5, 1, 5, 10, 20]

    # Create the Cask for the given reactor
    cask = Cask.from_reactor(reactor, total_mass=cask_mass)

    # Create the Spectra for the given removal times
    spectra = []
    for removal_time in removal_times:
        spec_05 = cask.get_total_spectrum(removal_time=removal_time)
        spectra.append(spec_05)

    # Compare to the reference data file
    # The reference data we have (saved by command_line.py) is for the
    # 0.5 year removal time only.
    spec_05 = spectra[1]
    with importlib.resources.path(data, f"{reactor.capitalize()}_single.csv") as path:
        spec_ref = Spectrum.from_file(path)

    assert len(spec_05.energy) == len(spec_ref.energy), (
        f"Energy data has wrong length "
        f"({len(spec_05.energy)} vs {len(spec_ref.energy)})"
    )
    assert len(spec_05.flux) == len(spec_ref.flux), (
        f"Flux data has wrong length ({len(spec_05.flux)} vs {len(spec_ref.flux)})"
    )
    assert np.allclose(spec_05.energy, spec_ref.energy), (
        "Reference energy does not match"
    )
    assert np.allclose(spec_05.flux, spec_ref.flux), "Reference flux does not match"

    # Calculate the single cask flux at 40m
    # Returns a single value of flux cm^-2 s^-1
    total_flux = spec_05.integrate(1806, 6000)
    flux_at_40m = calculate_flux_at_distance(total_flux, distance=40)
    assert isinstance(flux_at_40m, float), "Single flux is not a float"
    if reactor == "sizewell":
        flux_ref = 1191944053.3836975
    elif reactor == "hartlepool":
        flux_ref = 491107390.00662863
    assert np.isclose(flux_at_40m, flux_ref), (
        f"Single flux value does not match: ({flux_ref:e} vs {flux_at_40m:e})"
    )

    # Calculate event rates in a detector at 40m using the flux spectrum
    rate_lower, rate_upper = calculate_event_rate(flux_at_40m, 0.2, 0.4)
    assert isinstance(rate_lower, float), "Lower event rate is not a float"
    assert isinstance(rate_upper, float), "Upper event rate is not a float"
    if reactor == "sizewell":
        rate_lower_ref = 1.7734906963638752e-07
        rate_upper_ref = 3.5469813927277504e-07
    elif reactor == "hartlepool":
        rate_lower_ref = 7.307175069331267e-08
        rate_upper_ref = 1.4614350138662533e-07
    assert np.isclose(rate_lower, rate_lower_ref, atol=1e-15), (
        f"Lower event rate does not match: ({rate_lower_ref:e} vs {rate_lower:e})"
    )
    assert np.isclose(rate_upper, rate_upper_ref, atol=1e-15), (
        f"Upper event rate does not match: ({rate_upper_ref:e} vs {rate_upper:e})"
    )


@pytest.mark.parametrize("reactor", ["sizewell", "hartlepool"])
def test_multiple_casks(reactor: str) -> None:
    """Test multiple cask spectrum output."""
    cask_mass = 10000  # 10 tonne casks
    casks_per_removal = 10
    if reactor == "sizewell":
        removal_times = [0.5, 5, 10, 20]
    elif reactor == "hartlepool":
        removal_times = [3, 7, 15, 19]

    # Create the spectra and combine them
    cask = Cask.from_reactor(reactor, total_mass=cask_mass * casks_per_removal)
    spectra = []
    for removal_time in removal_times:
        spec = cask.get_total_spectrum(removal_time=removal_time)
        spectra.append(spec)
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
        flux_ref = 12982725680.40801
    elif reactor == "hartlepool":
        flux_ref = 1428527216.9547708
    assert np.isclose(flux_at_40m, flux_ref), (
        f"Total flux value does not match: ({flux_ref:e} vs {flux_at_40m:e})"
    )

    # Calculate event rates in a detector at 40m using the flux spectrum
    rate_lower, rate_upper = calculate_event_rate(flux_at_40m, 0.2, 0.4)
    assert isinstance(rate_lower, float), "Lower event rate is not a float"
    assert isinstance(rate_upper, float), "Upper event rate is not a float"
    if reactor == "sizewell":
        rate_lower_ref = 1.9316966381337443e-06
        rate_upper_ref = 3.863393276267489e-06
    elif reactor == "hartlepool":
        rate_lower_ref = 2.1255022176416828e-07
        rate_upper_ref = 4.2510044352833657e-07
    assert np.isclose(rate_lower, rate_lower_ref, atol=1e-15), (
        f"Lower event rate does not match: ({rate_lower_ref:e} vs {rate_lower:e})"
    )
    assert np.isclose(rate_upper, rate_upper_ref, atol=1e-15), (
        f"Upper event rate does not match: ({rate_upper_ref:e} vs {rate_upper:e})"
    )
