"""Basic functional tests to check the output of examples from the command line script.

These aren't 'real' tests yet, there area few issues like relative paths to data files
that need to be resolved before they can be fully integrated into a test suite.

For now this script has to be run from within the snf_simulations directory.

But at least it should check that the output from the command line script hasn't changed,
and give a baseline for future tests.
"""

import importlib.resources

import numpy as np
import pytest
import ROOT

from snf_simulations.physics import calculate_event_rate, calculate_flux
from snf_simulations.sample import sample_spec
from snf_simulations.scripts.command_line import _get_spectra
from snf_simulations.spec import add_spec

from . import test_data

ROOT.TH1.AddDirectory(False)  # Prevent ROOT from keeping histograms in memory

# Suppress warnings from ruff
# ruff: noqa: S101  # asserts


def _spec_to_arrays(spec):
    """Convert a ROOT.TH1D spectrum to energy and flux numpy arrays."""
    n_bins = spec.GetNbinsX()
    energy = np.array([spec.GetBinCenter(i) for i in range(1, n_bins + 1)])
    flux = np.array([spec.GetBinContent(i) for i in range(1, n_bins + 1)])
    return energy, flux


def _load_output(filename):
    """Load data from an output file.

    The 'CSV' files written by the old `flux.write_spec_single` and
    `flux.write_spec_multiple` functions are not really proper CSV files,
    so we have to do a bit of manual parsing here to get the data out.
    """
    with importlib.resources.path(test_data, filename) as path, open(path) as f:
        energy_line = f.readline()
        flux_line = f.readline()

    # Line format is:
    # "energy": [0.0005, 0.0015, ... 5.2995, 5.3005],\n
    # "(flux)": [58875835323384.055, ... 0.0, 0.0],\n
    energy = np.array([float(e) for e in energy_line[11:-3].split(", ")])
    flux = np.array([float(fx) for fx in flux_line[11:-3].split(", ")])

    # Energy was saved in MeV, convert to keV
    energy *= 1e3

    return energy, flux


def _load_output_new(filename):
    """Load data from a proper CSV output file."""
    with importlib.resources.path(test_data, filename) as path:
        data = np.loadtxt(
            path,
            delimiter=",",
            skiprows=1,
        )
    energy = data[:, 0].tolist()
    flux = data[:, 1].tolist()
    return energy, flux


@pytest.mark.parametrize("reactor", ["sizewell", "hartlepool"])
def test_single_cask(reactor):
    """Test single cask spectrum output."""
    cask_mass = 100000  # test data is for a 100 tonne cask (not 10 tonnes)
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
        assert len(counts) == len(values) == 5312, "Data has wrong length"

    # Compare to the reference data file
    # The reference data we have (saved by command_line.py) is for the
    # 0.5 year removal time only, and saved in a non-standard format.
    spec = spectra[1]
    energy_single, flux_single = _spec_to_arrays(spec)
    energy_ref, flux_ref = _load_output(f"{reactor.capitalize()}_single_0.5.csv")
    assert np.allclose(energy_ref, energy_single), "Reference energy does not match"
    assert np.allclose(flux_ref, flux_single), "Reference flux does not match"

    # Calculate the single cask flux at 40m
    # Returns a single value of flux cm^-2 s^-1
    flux = calculate_flux(spec, distance=40)
    assert isinstance(flux, float), "Single flux is not a float"
    if reactor == "sizewell":
        assert flux == 11992567783.00658, "Single flux value does not match"
    else:
        assert flux == 4942443091.633026, "Single flux value does not match"

    # Calculate event rates in a detector at 40m using the flux spectrum
    rate_lower, rate_upper = calculate_event_rate(flux, 0.2, 0.4)
    assert isinstance(rate_lower, float), "Lower event rate is not a float"
    assert isinstance(rate_upper, float), "Upper event rate is not a float"
    if reactor == "sizewell":
        assert rate_lower == 1.7843712822172811e-06, "Lower event rate does not match"
        assert rate_upper == 3.5687425644345623e-06, "Upper event rate does not match"
    else:
        assert rate_lower == 7.353849214177359e-07, "Lower event rate does not match"
        assert rate_upper == 1.4707698428354718e-06, "Upper event rate does not match"


@pytest.mark.parametrize("reactor", ["sizewell", "hartlepool"])
def test_multiple_casks(reactor):
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
    assert len(counts) == len(values) == 5312, "Data has wrong length"

    # Compare to the reference data file
    energy_multiple, flux_multiple = _spec_to_arrays(spec_multiple)
    energy_ref, flux_ref = _load_output(f"{reactor.capitalize()}_multiple.csv")
    assert np.allclose(energy_ref, energy_multiple), "Reference energy does not match"
    assert np.allclose(flux_ref, flux_multiple), "Reference flux does not match"

    # Calculate the total cask flux at 40m
    flux = calculate_flux(spec_multiple, distance=40)
    assert isinstance(flux, float), "Total flux is not a float"
    if reactor == "sizewell":
        assert flux == 13067897209.144945, "Total flux value does not match"
    else:
        assert flux == 1443268305.2325196, "Total flux value does not match"

    # Calculate event rates in a detector at 40m using the flux spectrum
    rate_lower, rate_upper = calculate_event_rate(flux, 0.2, 0.4)
    assert isinstance(rate_lower, float), "Lower event rate is not a float"
    assert isinstance(rate_upper, float), "Upper event rate is not a float"
    if reactor == "sizewell":
        assert rate_lower == 1.9443692894533463e-06, "Lower event rate does not match"
        assert rate_upper == 3.8887385789066926e-06, "Upper event rate does not match"
    else:
        assert rate_lower == 2.1474354475115334e-07, "Lower event rate does not match"
        assert rate_upper == 4.294870895023067e-07, "Upper event rate does not match"


@pytest.mark.parametrize("reactor", ["sizewell", "hartlepool"])
def test_sampling(reactor):
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
    samples = sample_spec(
        spec_multiple,
        samples=1000000,
    )

    # Check the samples
    assert len(samples) == 1000000, "Wrong number of samples in CSV"

    # Compare to reference file
    # sample.sample() uses GetRandom(), but the output seems to be deterministic.
    # So we can compare to a reference file.
    filename = f"{reactor.capitalize()}_sampled_spectrum.csv"
    with importlib.resources.path(test_data, filename) as path, open(path) as f:
        lines = f.readlines()
        samples_ref = [float(line.strip()) for line in lines]
    assert len(samples_ref) == 1000000, "Wrong number of samples in reference CSV"
    assert samples == samples_ref, "Sampled spectrum does not match reference"
