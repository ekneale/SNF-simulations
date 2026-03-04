"""Temporary tests between legacy ROOT functions and the new Numpy-based classes.

These tests will be skipped if ROOT is not installed, and will be removed
once the ROOT-based code is fully removed.
"""

import numpy as np
import pytest

try:
    import ROOT
except ImportError:
    pytest.skip("skipping ROOT tests", allow_module_level=True)

from snf_simulations.data import load_antineutrino_data, load_isotope_data
from snf_simulations.spec import Spectrum, equalise_spec

from .test_spec import _mock_data

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts

ROOT.TH1.AddDirectory(False)  # Prevent ROOT from keeping histograms in memory


_MOLAR_MASSES, _ = load_isotope_data()
ISOTOPES = list(_MOLAR_MASSES.keys())
ISOTOPE_SPECTRA = load_antineutrino_data(ISOTOPES)


def _get_mock_spectra() -> tuple[Spectrum, ROOT.TH1D]:
    """Create a mock Spectrum and corresponding ROOT TH1D histogram."""
    # Use the same mock data as in test_spec.
    energy, flux, errors = _mock_data()

    # Create the Spectrum instance
    # Note we miss off the final flux and error values, but still pass the full
    # array of energy values.
    spec = Spectrum(energy, flux[:-1], errors[:-1], name="mock")

    # Create the ROOT histogram
    # Note we exclude the final data row entirely.
    root_spec = ROOT.TH1D("mock", "mock", len(energy) - 1, np.array(energy))
    for e, f in zip(energy[:-1], flux[:-1], strict=True):
        root_spec.Fill(e, f)
    for i, err in enumerate(errors[:-1], start=1):
        root_spec.SetBinError(i, err)

    return spec, root_spec


def _root_to_arrays(spec: ROOT.TH1D) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert a ROOT TH1D histogram to Numpy arrays."""
    n_bins = spec.GetNbinsX()
    lower_edges = [spec.GetBinLowEdge(i) for i in range(1, n_bins + 1)]
    upper_edge = spec.GetBinLowEdge(n_bins) + spec.GetBinWidth(n_bins)
    edges = np.array([*lower_edges, upper_edge], dtype=float)
    content = np.array(
        [spec.GetBinContent(i) for i in range(1, n_bins + 1)],
        dtype=float,
    )
    errors = np.array([spec.GetBinError(i) for i in range(1, n_bins + 1)], dtype=float)
    return edges, content, errors


def _compare_spectra(spec: Spectrum, root_spec: ROOT.TH1D) -> None:
    """Compare a Spectrum instance to a ROOT TH1D histogram."""
    root_energy, root_flux, root_errors = _root_to_arrays(root_spec)

    assert len(spec.energy) == len(root_energy), (
        f"Spectrum '{spec.name}' energy array has different length to ROOT energy array"
        f" ({len(spec.energy)} vs {len(root_energy)})"
    )
    assert len(spec.flux) == len(root_flux), (
        f"Spectrum '{spec.name}' flux array has different length to ROOT flux array"
        f" ({len(spec.flux)} vs {len(root_flux)})"
    )
    assert len(spec.errors) == len(root_errors), (
        f"Spectrum '{spec.name}' error array has different length to ROOT error array"
        f" ({len(spec.errors)} vs {len(root_errors)})"
    )

    assert np.allclose(spec.energy, root_energy), (
        f"Spectrum '{spec.name}' energy edges do not match ROOT edges"
    )
    assert np.allclose(spec.flux, root_flux), (
        f"Spectrum '{spec.name}' flux values do not match ROOT contents"
    )
    assert np.allclose(spec.errors, root_errors), (
        f"Spectrum '{spec.name}' error values do not match ROOT errors"
    )


def test_spectrum() -> None:
    """Test that the Spectrum class matches ROOT."""
    spec, root_spec = _get_mock_spectra()

    # Compare the arrays to ensure they match
    _compare_spectra(spec, root_spec)


def test_equalise() -> None:
    """Test that Spectrum.equalise matches ROOT."""
    spec, root_spec = _get_mock_spectra()

    # Equalise the spectra
    spec.equalise(width=1, min_energy=0, max_energy=3)
    root_equal = equalise_spec(root_spec, max_energy=3, min_energy=0)

    # Compare the arrays to ensure they match
    _compare_spectra(spec, root_equal)


def test_scale() -> None:
    """Test that Spectrum scaling matches ROOT."""
    spec, root_spec = _get_mock_spectra()

    # Scale the spectra
    scale_factor = 2.5
    spec_scaled = spec * scale_factor
    root_scaled = root_spec.Clone()
    root_scaled.Scale(scale_factor)

    # Compare the arrays to ensure they match
    _compare_spectra(spec_scaled, root_scaled)


def test_add() -> None:
    """Test that adding two Spectra matches adding two ROOT spectra."""
    spec1, root_spec1 = _get_mock_spectra()
    spec2, root_spec2 = _get_mock_spectra()

    # Add the spectra together
    spec_sum = spec1 + spec2
    root_spectra = ROOT.TList()
    root_spectra.Add(root_spec1)
    root_spectra.Add(root_spec2)
    root_sum = root_spectra[0].Clone("combined")
    root_sum.Reset()
    root_sum.Merge(root_spectra)

    # Compare the arrays to ensure they match
    _compare_spectra(spec_sum, root_sum)

    # Check that the sum of two identical spectra is double the original
    spec_scaled = spec1 * 2
    root_scaled = root_spec1.Clone()
    root_scaled.Scale(2)
    _compare_spectra(spec_scaled, root_scaled)


def test_sample() -> None:
    """Test that Spectrum.sample produces a distribution comparable to ROOT's."""
    spec, root_spec = _get_mock_spectra()

    # Sample from both spectra
    samples = 100000
    spec_samples = spec.sample(samples=samples)
    root_samples = np.array([root_spec.GetRandom() for _ in range(samples)])
    assert len(spec_samples) == samples, (
        "Spectrum.sample did not produce the expected number of samples"
    )
    assert len(root_samples) == samples, (
        "ROOT.GetRandom did not produce the expected number of samples"
    )

    # We can't compare the samples directly due to the random selection, but we can
    # compare their distributions to the expected probabilities from the histogram.
    # First, compare against the expected probabilities which are proportional to
    # the (normalised) flux values in each bin.
    expected_probabilities = spec.flux / np.sum(spec.flux)
    spec_counts, _ = np.histogram(spec_samples, bins=spec.energy)
    root_counts, _ = np.histogram(root_samples, bins=spec.energy)

    # We allow a reasonable 5-sigma error margin, where sigma is the standard deviation
    # of the binomial distribution for each bin.
    sigma = np.sqrt(expected_probabilities * (1 - expected_probabilities) / samples)
    spec_probabilities = spec_counts / samples
    root_probabilities = root_counts / samples
    assert np.all(np.abs(spec_probabilities - expected_probabilities) < 5 * sigma), (
        "Spectrum.sample bin probabilities deviate more than expected from "
        "the target histogram distribution"
    )
    assert np.all(np.abs(root_probabilities - expected_probabilities) < 5 * sigma), (
        "ROOT.GetRandom bin probabilities deviate more than expected from "
        "the target histogram distribution"
    )

    # Two-sample chi-squared consistency check between the two outputs.
    combined_counts = root_counts + spec_counts
    # Exclude bins where both counts are zero to avoid division by zero in chi2.
    valid = combined_counts > 0
    chi2 = np.sum(
        (root_counts[valid] - spec_counts[valid]) ** 2 / combined_counts[valid],
    )
    degrees_of_freedom = int(np.count_nonzero(valid) - 1)
    assert degrees_of_freedom > 0

    # The standard deviation of the chi-squared distribution is sqrt(2*ndf),
    # again allow a 5-sigma deviation from the mean.
    sigma = np.sqrt(2 * degrees_of_freedom)
    assert chi2 < degrees_of_freedom + 5 * sigma, (
        f"ROOT and Spectrum sampled distributions differ too much: "
        f"chi2={chi2:.3f}, ndf={degrees_of_freedom}"
    )


def test_integrate() -> None:
    """Test that Spectrum.integrate works for arbitrary (non-equalised) bins."""
    spec, root_spec = _get_mock_spectra()

    # Test integration over full range.
    # Flux is dN/dE (keV^-1), so integral is sum(flux * bin_width).
    expected_integral = np.sum(spec.flux * np.diff(spec.energy))
    spec_integral = spec.integrate(spec.energy[0], spec.energy[-1])
    root_integral = root_spec.Integral(1, root_spec.GetNbinsX(), "width")

    assert np.isclose(spec_integral, expected_integral), (
        f"Spectrum.integrate({spec.energy[0]}, {spec.energy[-1]}) = {spec_integral}"
        f" does not match expected value {expected_integral}"
    )
    assert np.isclose(root_integral, expected_integral), (
        f"ROOT spec.Integral(1, {root_spec.GetNbinsX()}, 'width') = {root_integral}"
        f" does not match expected value {expected_integral}"
    )

    # Test integration over partial range (still full bins).
    lower_energy = 0.5
    upper_energy = 3.0
    lower_index = 1
    upper_index = 4
    expected_integral = np.sum(
        spec.flux[lower_index:upper_index]
        * np.diff(spec.energy)[lower_index:upper_index],
    )
    spec_integral = spec.integrate(lower_energy, upper_energy)
    root_integral = root_spec.Integral(lower_index + 1, upper_index, "width")
    assert np.isclose(spec_integral, expected_integral), (
        f"Spectrum.integrate({lower_energy}, {upper_energy}) = {spec_integral}"
        f" does not match expected value {expected_integral}"
    )
    assert np.isclose(root_integral, expected_integral), (
        f"ROOT spec.Integral({lower_index + 1}, {upper_index}, 'width') = "
        f"{root_integral}"
        f" does not match expected value {expected_integral}"
    )

    # Currently, we only ever use ROOT's Integral on spectra that have been equalised
    # to 1 keV bins from 0 to 6000. That means we don't have to worry about partial
    # bin widths and can use the energy range values as bin values.
    # For the Spectrum class we want to be able to integrate over arbitrary energy
    # ranges, including partial bins.
    # We can test that elsewhere, since trying to get ROOT to do partial bin integration
    # is non-trivial.


@pytest.mark.parametrize("isotope", ISOTOPES)
def test_spec_real(isotope: str) -> None:
    """Test that Spectra from real data matches ROOT."""
    max_energy = 6000
    isotope_data = ISOTOPE_SPECTRA[isotope]
    energy = isotope_data[:, 0]
    flux = isotope_data[:, 1]
    errors = isotope_data[:, 2]

    # Create the ROOT spectrum and equalise it
    root_spec = ROOT.TH1D(isotope, isotope, len(energy) - 1, np.array(energy))
    for e, f in zip(energy[:-1], flux[:-1], strict=True):
        root_spec.Fill(e, f)
    for i, err in enumerate(errors[:-1], start=1):
        root_spec.SetBinError(i, err)
    root_equal = equalise_spec(root_spec, max_energy=max_energy, min_energy=0)

    # Create the Spectrum instance and equalise it
    spec = Spectrum.from_isotope(isotope)
    spec.equalise(width=1, min_energy=0, max_energy=max_energy)

    # Compare the arrays to ensure they match
    _compare_spectra(spec, root_equal)
