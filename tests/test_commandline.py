"""Unit tests for command line script functions."""

import importlib.resources

import numpy as np
import pytest
import ROOT

from snf_simulations.physics import calculate_event_rate, calculate_flux_at_distance
from snf_simulations.scripts.command_line import _get_spectra
from snf_simulations.spec import add_spec, sample_spec

from . import data

ROOT.TH1.AddDirectory(False)  # Prevent ROOT from keeping histograms in memory

# Suppress warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers


def _spec_to_arrays(spec: ROOT.TH1D) -> tuple[np.ndarray, np.ndarray]:
    """Convert a ROOT.TH1D spectrum to energy and flux numpy arrays."""
    n_bins = spec.GetNbinsX()
    energy = np.array([spec.GetBinCenter(i) for i in range(1, n_bins + 1)])
    flux = np.array([spec.GetBinContent(i) for i in range(1, n_bins + 1)])
    return energy, flux


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

    # Create the spectra for the given removal times
    spectra = _get_spectra(
        cask_mass=cask_mass,
        removal_times=removal_times,
        reactor=reactor,
    )

    # Check the returned list of spectra
    assert isinstance(spectra, ROOT.TList), "Returned object is not a TList"
    assert len(spectra) == len(removal_times), "Returned list has wrong length"
    assert all(isinstance(spec, ROOT.TH1D) for spec in spectra), (
        "Not all items in list are TH1D"
    )

    # Check the basic properties of each of the spectra
    for spec in spectra:
        assert isinstance(spec, ROOT.TH1D), "Returned spectrum is not a TH1D"
        counts, values = spec.counts(), spec.values()
        assert isinstance(counts, np.ndarray), "Counts are not a numpy array"
        assert isinstance(values, np.ndarray), "Values are not a numpy array"

    # Compare to the reference data file
    # The reference data we have (saved by command_line.py) is for the
    # 0.5 year removal time only.
    spec = spectra[1]
    energy_single, flux_single = _spec_to_arrays(spec)
    energy_ref, flux_ref = _load_output(f"{reactor.capitalize()}_single.csv")

    assert len(energy_ref) == len(energy_single), (
        f"Energy data has wrong length ({len(energy_ref)} vs {len(energy_single)})"
    )
    assert len(flux_ref) == len(flux_single), (
        f"Flux data has wrong length ({len(flux_ref)} vs {len(flux_single)})"
    )
    assert np.allclose(energy_ref, energy_single), "Reference energy does not match"
    assert np.allclose(flux_ref, flux_single), "Reference flux does not match"

    # Calculate the single cask flux at 40m
    # Returns a single value of flux cm^-2 s^-1
    total_flux = spec.Integral(1806, 6000)
    flux_at_40m = calculate_flux_at_distance(total_flux, distance=40)
    assert isinstance(flux_at_40m, float), "Single flux is not a float"
    if reactor == "sizewell":
        flux_ref = 1193212153.241669
    elif reactor == "hartlepool":
        flux_ref = 491654259.98276895
    assert np.isclose(flux_at_40m, flux_ref), (
        f"Single flux value does not match: ({flux_ref:e} vs {flux_at_40m:e})"
    )

    # Calculate event rates in a detector at 40m using the flux spectrum
    rate_lower, rate_upper = calculate_event_rate(flux_at_40m, 0.2, 0.4)
    assert isinstance(rate_lower, float), "Lower event rate is not a float"
    assert isinstance(rate_upper, float), "Upper event rate is not a float"
    if reactor == "sizewell":
        rate_lower_ref = 1.7848099533254326e-07
        rate_upper_ref = 3.569619906650865e-07
    elif reactor == "hartlepool":
        rate_lower_ref = 7.355992829950914e-08
        rate_upper_ref = 1.4711985659901828e-07
    assert np.isclose(rate_lower, rate_lower_ref), (
        f"Lower event rate does not match: ({rate_lower_ref:e} vs {rate_lower:e})"
    )
    assert np.isclose(rate_upper, rate_upper_ref), (
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
    spectra = _get_spectra(
        cask_mass=cask_mass * casks_per_removal,
        removal_times=removal_times,
        reactor=reactor,
    )
    spec_multiple = add_spec(spectra)

    # Check the combined spectrum
    assert isinstance(spec_multiple, ROOT.TH1D), "Combined spectrum is not a TH1D"
    counts, values = spec_multiple.counts(), spec_multiple.values()
    assert isinstance(counts, np.ndarray), "Counts are not a numpy array"
    assert isinstance(values, np.ndarray), "Values are not a numpy array"

    # Compare to the reference data file
    energy_multiple, flux_multiple = _spec_to_arrays(spec_multiple)
    energy_ref, flux_ref = _load_output(f"{reactor.capitalize()}_multiple.csv")
    assert len(energy_ref) == len(energy_multiple), (
        f"Energy data has wrong length ({len(energy_ref)} vs {len(energy_multiple)})"
    )
    assert len(flux_ref) == len(flux_multiple), (
        f"Flux data has wrong length ({len(flux_ref)} vs {len(flux_multiple)})"
    )
    assert np.allclose(energy_ref, energy_multiple), "Reference energy does not match"
    assert np.allclose(flux_ref, flux_multiple), "Reference flux does not match"

    # Calculate the total cask flux at 40m
    total_flux = spec_multiple.Integral(1806, 6000)
    flux_at_40m = calculate_flux_at_distance(total_flux, distance=40)
    assert isinstance(flux_at_40m, float), "Total flux is not a float"
    if reactor == "sizewell":
        flux_ref = 12997598681.459194
    elif reactor == "hartlepool":
        flux_ref = 1431196201.737153
    assert np.isclose(flux_at_40m, flux_ref), (
        f"Total flux value does not match: ({flux_ref:e} vs {flux_at_40m:e})"
    )

    # Calculate event rates in a detector at 40m using the flux spectrum
    rate_lower, rate_upper = calculate_event_rate(flux_at_40m, 0.2, 0.4)
    assert isinstance(rate_lower, float), "Lower event rate is not a float"
    assert isinstance(rate_upper, float), "Upper event rate is not a float"
    if reactor == "sizewell":
        rate_lower_ref = 1.93390958839063e-06
        rate_upper_ref = 3.86781917678126e-06
    elif reactor == "hartlepool":
        rate_lower_ref = 2.1294733936938256e-07
        rate_upper_ref = 4.2589467873876513e-07
    assert np.isclose(rate_lower, rate_lower_ref), (
        f"Lower event rate does not match: ({rate_lower_ref:e} vs {rate_lower:e})"
    )
    assert np.isclose(rate_upper, rate_upper_ref), (
        f"Upper event rate does not match: ({rate_upper_ref:e} vs {rate_upper:e})"
    )


@pytest.mark.parametrize("reactor", ["sizewell", "hartlepool"])
def test_sampling(reactor: str) -> None:
    """Test plotting a sample spectrum for given removal times."""
    cask_mass = 10000  # 10 tonne casks
    casks_per_removal = 10
    if reactor == "sizewell":
        removal_times = [0.5, 5, 10, 20]
    elif reactor == "hartlepool":
        removal_times = [3, 7, 15, 19]

    # Create the spectra and combine them
    spectra = _get_spectra(
        cask_mass=cask_mass * casks_per_removal,
        removal_times=removal_times,
        reactor=reactor,
    )
    spec_multiple = add_spec(spectra)

    # Sample the spectrum
    n_samples = 1000000
    samples = sample_spec(
        spec_multiple,
        samples=n_samples,
    )

    # Check the samples
    assert len(samples) == n_samples, "Wrong number of samples selected"

    # We can't compare the samples directly due to the random selection, but we can
    # compare their distributions to the expected probabilities from the histogram.
    # First, compare against the expected probabilities which are proportional to
    # the (normalised) flux values in each bin.
    _, flux = _spec_to_arrays(spec_multiple)
    expected_probabilities = flux / np.sum(flux)
    n_bins = spec_multiple.GetNbinsX()
    lower_edges = [spec_multiple.GetBinLowEdge(i) for i in range(1, n_bins + 1)]
    upper_edge = spec_multiple.GetBinLowEdge(n_bins) + spec_multiple.GetBinWidth(n_bins)
    bin_edges = np.array([*lower_edges, upper_edge], dtype=float)
    counts, _ = np.histogram(samples, bins=bin_edges)

    # Use a chi-squared goodness-of-fit check against the expected bin probabilities.
    expected_counts = expected_probabilities * n_samples
    # Exclude bins where counts are zero to avoid division by zero in chi2.
    valid = expected_counts > 0
    chi2 = np.sum(
        (counts[valid] - expected_counts[valid]) ** 2 / expected_counts[valid]
    )
    degrees_of_freedom = int(np.count_nonzero(valid) - 1)
    assert degrees_of_freedom > 0

    # The standard deviation of the chi-squared distribution is sqrt(2*ndf),
    # again allow a 5-sigma deviation from the mean.
    sigma = np.sqrt(2 * degrees_of_freedom)
    assert chi2 < degrees_of_freedom + 5 * sigma, (
        f"Sampled distribution differs too much from expected histogram: "
        f"chi2={chi2:.3f}, ndf={degrees_of_freedom}"
    )
