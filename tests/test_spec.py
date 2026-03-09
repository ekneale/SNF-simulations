"""Unit tests for loading antineutrino spectra data."""

from pathlib import Path

import numpy as np
import pytest

from snf_simulations.data import load_antineutrino_data, load_isotope_data
from snf_simulations.spec import Spectrum

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers

_MOLAR_MASSES, _ = load_isotope_data()
ISOTOPES = list(_MOLAR_MASSES.keys())


def _mock_data() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a standard mock spectrum dataset used across tests."""
    # Example of energy, flux (dN/dE), and uncertainty data for a mock spectrum.
    # Note the deliberately uneven energy bin spacing.
    data = np.array(
        [
            [0.0, 10, 1],
            [0.5, 20, 2],
            [1.0, 30, 3],
            [2.0, 40, 4],
            [3.0, 50, 5],
            [5.0, 60, 6],
            [8.0, 70, 7],
        ],
    )
    energy = data[:, 0]
    flux = data[:, 1]
    errors = data[:, 2]
    return energy, flux, errors


def test_load_isotope_data() -> None:
    """Test that isotope data is imported correctly."""
    molar_masses, half_lives = load_isotope_data()
    assert isinstance(molar_masses, dict), "Molar masses is not a dictionary"
    assert isinstance(half_lives, dict), "Half lives is not a dictionary"
    assert len(molar_masses) == 16, f"Loaded {len(molar_masses)} isotopes, expected 16"
    assert len(molar_masses) == len(half_lives), (
        "Molar masses and half lives have different number of isotopes",
    )


def test_load_antineutrino_data() -> None:
    """Test that antineutrino spectra are imported correctly."""
    molar_masses, _ = load_isotope_data()
    isotopes = list(molar_masses.keys())
    data = load_antineutrino_data(isotopes)
    assert isinstance(data, dict), "Loaded data is not a dictionary"
    assert len(data) == 16, f"Loaded {len(data)} isotopes, expected 16"

    for isotope, isotope_data in data.items():
        # Basic structure checks
        assert isinstance(isotope_data, np.ndarray), (
            f"Isotope {isotope} spectrum is not a numpy array"
        )
        assert isotope_data.shape[1] == 3, (
            f"Isotope {isotope} spectrum has {isotope_data.shape[1]} columns,"
            f" expected 3"
        )
        assert isotope_data.shape[0] > 0, f"Isotope {isotope} spectrum has no rows"

        # Physical value checks
        assert np.all(isotope_data[:, 0] >= 0), (
            f"Isotope {isotope} has negative energies"
        )
        assert np.all(isotope_data[:, 1] >= 0), (
            f"Isotope {isotope} has negative flux values"
        )
        assert np.all(isotope_data[:, 2] >= 0), (
            f"Isotope {isotope} has negative uncertainties"
        )


def test_create_spec() -> None:
    """Test Spectrum construction with mock data."""
    # Create the Spectrum
    # Note we remove the last flux & error values, since the last energy value is
    # taken as the upper edge of the penultimate bin, and we don't have an upper edge
    # for the final flux bin.
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    # Test basic properties
    assert isinstance(spec, Spectrum), "Spectrum is not a Spectrum"
    assert len(spec.energy) == len(energy), (
        f"Spectrum has {len(spec.energy)} energy values instead of {len(energy)}"
    )
    assert len(spec.flux) == len(energy) - 1, (
        f"Spectrum has {len(spec.flux)} flux values instead of {len(energy) - 1}"
    )
    assert len(spec.errors) == len(energy) - 1, (
        f"Spectrum has {len(spec.errors)} error values instead of {len(energy) - 1}"
    )
    assert np.array_equal(spec.energy, energy), "Spectrum energy values mismatch"
    assert np.array_equal(spec.flux, flux[:-1]), "Spectrum flux values mismatch"
    assert np.array_equal(spec.errors, errors[:-1]), "Spectrum error values mismatch"


def test_init() -> None:
    """Test that Spectrum constructor validates input shapes and lengths."""
    energy, flux, errors = _mock_data()

    with pytest.raises(ValueError, match="Energy, flux and errors must be 1D arrays"):
        Spectrum(energy.reshape(-1, 1), flux, errors)

    with pytest.raises(
        ValueError,
        match=r"Flux array must have length len\(energy\) - 1",
    ):
        Spectrum(energy, flux, errors)

    with pytest.raises(
        ValueError,
        match="Flux and errors arrays must have the same length",
    ):
        Spectrum(energy, flux[:-1], errors[:-2])


def test_repr() -> None:
    """Test repr for initialized and uninitialized Spectrum objects."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1], name="mock")

    expected_repr = (
        f'<Spectrum "{spec.name}", '
        f"energy_range=({spec.energy[0]:.1f}-{spec.energy[-1]:.1f} keV)>"
    )
    assert repr(spec) == expected_repr, (
        "Spectrum repr should include spectrum name and energy range"
    )

    # Create an instance without calling __init__ to exercise fallback repr path.
    uninitialized = Spectrum.__new__(Spectrum)
    assert repr(uninitialized) == "<Spectrum (uninitialized)>", (
        "Spectrum repr should handle uninitialized objects gracefully"
    )


@pytest.mark.parametrize("isotope", ISOTOPES)
def test_from_isotope(isotope: str) -> None:
    """Test that Spectrum.from_isotope returns a valid spectrum object."""
    spec = Spectrum.from_isotope(isotope)
    assert isinstance(spec, Spectrum), "from_isotope should return a Spectrum"
    assert spec.name == isotope, "from_isotope should preserve isotope name"
    assert len(spec.energy) == len(spec.flux) + 1, "Energy edges/flux length mismatch"
    assert np.all(np.diff(spec.energy) > 0), (
        "Energy edges should be strictly increasing"
    )


@pytest.mark.parametrize("width", [0.1, 1.0, 2.0])
def test_equalise(width: float) -> None:
    """Test Spectrum equalisation over the energy range for given bin width."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    # Equalise the Spectrum between the min and max energy values
    min_energy = int(np.min(energy))
    max_energy = int(np.max(energy))
    spec.equalise(width=width, min_energy=min_energy, max_energy=max_energy)

    # Test basic properties
    expected_bins = int((max_energy - min_energy) / width)
    assert len(spec.energy) == expected_bins + 1, (
        f"Spectrum has {len(spec.energy)} energy values instead of {expected_bins + 1}"
    )
    assert len(spec.flux) == expected_bins, (
        f"Spectrum has {len(spec.flux)} flux values instead of {expected_bins}"
    )
    assert len(spec.errors) == expected_bins, (
        f"Spectrum has {len(spec.errors)} error values instead of {expected_bins}"
    )

    # Test bin spacing
    # The bins should now be equally spaced integers from min_energy to max_energy.
    expected_edges = np.arange(min_energy, max_energy + width, width)
    assert np.all(spec.energy == expected_edges), (
        "Spectrum energy values mismatch after equalisation"
    )


@pytest.mark.parametrize("width", [0.1, 1.0, 2.0])
def test_equalise_extrapolate(width: float) -> None:
    """Test equalisation when extrapolating beyond original range."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    # Equalise the Spectrum outside of the current range
    min_energy = int(np.min(energy)) - 2
    max_energy = int(np.max(energy)) + 2
    spec.equalise(width=width, min_energy=min_energy, max_energy=max_energy)

    # Test bin spacing
    edges = np.arange(min_energy, max_energy + width, width)
    assert np.all(spec.energy == edges), (
        "Spectrum bin edges mismatch after equalisation with extrapolation"
    )

    # Test extrapolated values
    # Any bins fully outside the original range should be zero
    lower_edges = spec.energy[:-1]
    upper_edges = spec.energy[1:]
    outside = (upper_edges <= energy[0]) | (lower_edges >= energy[-1])
    assert np.all(spec.flux[outside] == 0), "Extrapolated bins should have zero flux"
    assert np.all(spec.errors[outside] == 0), (
        "Extrapolated bins should have zero errors"
    )


def test_equalise_extrapolate_overlap() -> None:
    """Test extrapolated bins that overlap with the original energy range."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    # This creates one fully outside low-energy bin [-0.2, 0.0]
    # and one overlapping bin [0.0, 0.2] that should remain non-zero.
    spec.equalise(width=0.2, min_energy=-0.2, max_energy=0.6)

    lower_edges = spec.energy[:-1]
    upper_edges = spec.energy[1:]
    outside = (upper_edges <= energy[0]) | (lower_edges >= energy[-1])
    overlap = ~outside
    assert np.all(spec.flux[outside] == 0), "Extrapolated bins should have zero flux"
    assert np.all(spec.errors[outside] == 0), (
        "Extrapolated bins should have zero errors"
    )
    assert np.any(spec.flux[overlap] > 0), (
        "Overlapping bins should retain non-zero flux"
    )


def test_equalise_integral() -> None:
    """Test that extrapolating does not change total integral."""
    energy, flux, errors = _mock_data()

    # Equalise the Spectrum between the min and max energy values
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])
    min_energy = int(np.min(energy))
    max_energy = int(np.max(energy))
    spec.equalise(width=1, min_energy=min_energy, max_energy=max_energy)

    # Equalise the Spectrum outside of the current range
    spec_extrapolated = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])
    spec_extrapolated.equalise(
        width=1, min_energy=min_energy - 2, max_energy=max_energy + 2
    )

    # Check that the integrals are the same
    base_integral = spec.integrate(spec.energy[0], spec.energy[-1])
    extrapolated_integral = spec_extrapolated.integrate(
        spec_extrapolated.energy[0], spec_extrapolated.energy[-1]
    )
    assert np.isclose(extrapolated_integral, base_integral), (
        "Extrapolating equalised range should not change total integral"
    )


def test_equalise_inputs() -> None:
    """Test equalise validation and default min/max energy behavior."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    with pytest.raises(ValueError, match="width must be a positive value"):
        spec.equalise(width=0)

    with pytest.raises(ValueError, match="max_energy must be greater than min_energy"):
        spec.equalise(width=1, min_energy=2, max_energy=2)

    spec.equalise(width=1)
    assert spec.energy[0] == energy[0], (
        "Min energy should default to original min energy when not specified"
    )
    assert spec.energy[-1] == energy[-1], (
        "Max energy should default to original max energy when not specified"
    )

    with pytest.raises(
        ValueError, match="width is too large for the given energy range"
    ):
        spec.equalise(width=10, min_energy=0, max_energy=8)


def test_scale() -> None:
    """Test positive scaling of flux and errors."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    # Scale the Spectrum
    scale_factor = 2.0
    scaled_spec = spec * scale_factor

    # Test basic properties
    assert isinstance(scaled_spec, Spectrum), "Scaled spectrum is not a Spectrum"
    assert np.all(np.isclose(scaled_spec.energy, spec.energy)), (
        "Scaled spectrum has different energy values than original spectrum"
    )

    # Test that each bin content and error has been scaled correctly
    expected_flux = spec.flux * scale_factor
    expected_errors = spec.errors * scale_factor
    assert np.all(np.isclose(scaled_spec.flux, expected_flux)), (
        "Scaled spectrum flux values mismatch expected values"
    )
    assert np.all(np.isclose(scaled_spec.errors, expected_errors)), (
        "Scaled spectrum error values mismatch expected values"
    )


def test_scale_negative() -> None:
    """Test scaling by a negative factor flips flux but keeps errors positive."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    scaled = spec * -3
    assert np.array_equal(scaled.flux, -3 * spec.flux), (
        "Negative scaling should flip flux sign"
    )
    assert np.array_equal(scaled.errors, np.abs(-3 * spec.errors)), (
        "Errors should scale with absolute factor"
    )


def test_add() -> None:
    """Test adding two spectra with matching binning."""
    energy, flux, errors = _mock_data()
    spec1 = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])
    spec2 = Spectrum(energy=energy, flux=flux[:-1] * 2, errors=errors[:-1] * 2)

    # Add the spectra together
    combined_spec = spec1 + spec2

    # Test basic properties
    assert isinstance(combined_spec, Spectrum), "Combined spectrum is not a Spectrum"
    assert len(combined_spec.energy) == len(spec1.energy), (
        "Combined spectrum has different number of energy bins"
    )

    # Test that each bin content and error has been added correctly
    expected_flux = spec1.flux + spec2.flux
    expected_errors = np.sqrt(spec1.errors**2 + spec2.errors**2)
    assert np.all(np.isclose(combined_spec.flux, expected_flux)), (
        "Combined spectrum flux values mismatch expected values"
    )
    assert np.all(np.isclose(combined_spec.errors, expected_errors)), (
        "Combined spectrum error values mismatch expected values"
    )


def test_add_different_energy() -> None:
    """Test that adding spectra with different binning raises an error."""
    energy, flux, errors = _mock_data()
    spec1 = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])
    spec2 = Spectrum(energy=energy * 2, flux=flux[:-1] * 2, errors=errors[:-1] * 2)

    with pytest.raises(
        ValueError,
        match="Energy bins of the two spectra must be the same",
    ):
        _ = spec1 + spec2


def test_sample_distribution() -> None:
    """Test deterministic sampling for a fixed seed."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])
    n_samples = 100000

    # Sample the spectrum
    samples = spec.sample(samples=n_samples)
    assert len(samples) == n_samples, "Wrong number of samples selected"

    # We can't compare the samples directly due to the random selection, but we can
    # compare their distributions to the expected probabilities from the histogram.
    # First, compare against the expected probabilities which are proportional to
    # the (normalised) flux values in each bin.
    expected_probabilities = spec.flux / np.sum(spec.flux)
    counts, _ = np.histogram(samples, bins=spec.energy)

    # Use a chi-squared goodness-of-fit check against the expected bin probabilities.
    expected_counts = expected_probabilities * n_samples
    # Exclude bins where counts are zero to avoid division by zero in chi2.
    valid = expected_counts > 0
    chi2 = np.sum(
        (counts[valid] - expected_counts[valid]) ** 2 / expected_counts[valid]
    )
    degrees_of_freedom = int(np.count_nonzero(valid) - 1)
    assert degrees_of_freedom > 0, (
        "Degrees of freedom must be positive for chi-squared test"
    )

    # The standard deviation of the chi-squared distribution is sqrt(2*ndf),
    # again allow a 5-sigma deviation from the mean.
    sigma = np.sqrt(2 * degrees_of_freedom)
    assert chi2 < degrees_of_freedom + 5 * sigma, (
        f"Sampled distribution differs too much from expected histogram: "
        f"chi2={chi2:.3f}, ndf={degrees_of_freedom}"
    )


def test_sample_with_seed() -> None:
    """Test deterministic sampling for a fixed seed."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    # Test sampling 5 values from the spectrum
    n_samples = 5
    samples = spec.sample(samples=n_samples, seed=1234)  # Fix seed for reproducibility

    # Samples should be fixed given the seed
    samples_ref = [
        5.3542737,
        2.24176629,
        5.95560179,
        1.96407925,
        2.2636498,
    ]
    assert len(samples) == n_samples, (
        f"Number of samples does not match requested number"
        f" ({len(samples)} vs {n_samples})"
    )
    assert np.allclose(samples, samples_ref), "Sampled values do not match reference"


def test_sample_reproducibility() -> None:
    """Test sampling reproducibility with constant seed."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    # Test reproducibility with a random (but constant) seed
    rng = np.random.default_rng()
    seed = int(rng.integers(0, int(1e6)))
    samples1 = spec.sample(samples=2000, seed=seed)
    samples2 = spec.sample(samples=2000, seed=seed)
    assert np.array_equal(samples1, samples2), (
        "Sampling should be reproducible with constant seed"
    )


def test_sample_zero_flux() -> None:
    """Test that sampling raises an error when total histogram weight is zero."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=np.zeros(len(flux[:-1])), errors=errors[:-1])

    with pytest.raises(
        ValueError,
        match="Histogram has zero total area; cannot sample",
    ):
        _ = spec.sample(samples=10)


def test_integrate() -> None:
    """Test integration over full and partial energy ranges."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    # Test integration over the full range
    min_energy = energy[0]
    max_energy = energy[-1]
    total_integral = spec.integrate(min_energy, max_energy)
    expected_total = np.sum(flux[:-1] * np.diff(energy))
    assert np.isclose(total_integral, expected_total), (
        "Total integral should match sum of bin content times bin width"
    )

    # Test integration over a partial range with full bins
    min_energy = energy[1]
    max_energy = energy[4]
    partial_integral = spec.integrate(min_energy, max_energy)
    expected_partial = np.sum(flux[:-1][1:4] * np.diff(energy)[1:4])
    assert np.isclose(partial_integral, expected_partial), (
        "Partial integral should include fractional contributions from edge bins"
    )

    # Test integration over a partial range that includes fractional bins
    min_energy = energy[1] - (energy[1] - energy[0]) * 0.5
    max_energy = energy[-2] + (energy[-1] - energy[-2]) * 0.5
    partial_integral = spec.integrate(min_energy, max_energy)
    expected_partial = (
        0.5 * flux[0] * (energy[1] - energy[0])  # half of the first bin
        + np.sum(flux[:-1][1:-1] * np.diff(energy)[1:-1])  # full middle bins
        + 0.5 * flux[-2] * (energy[-1] - energy[-2])  # half of the last bin
    )
    assert np.isclose(partial_integral, expected_partial), (
        "Partial integral should include fractional contributions from edge bins"
    )


def test_integrate_defaults() -> None:
    """Test that integrate defaults to full energy range."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    # Test both defaults together
    total_integral = spec.integrate()
    expected_total = spec.integrate(energy[0], energy[-1])
    assert np.isclose(total_integral, expected_total), (
        "Total integral should match sum of bin content times bin width"
    )

    # Test lower energy default
    partial_integral = spec.integrate(upper_energy=energy[-2])
    expected_partial = spec.integrate(energy[0], energy[-2])
    assert np.isclose(partial_integral, expected_partial), (
        "Partial integral with default lower energy should match expected value"
    )

    # Test upper energy default
    partial_integral = spec.integrate(lower_energy=energy[1])
    expected_partial = spec.integrate(energy[1], energy[-1])
    assert np.isclose(partial_integral, expected_partial), (
        "Partial integral with default upper energy should match expected value"
    )


def test_integrate_outside_range() -> None:
    """Test that integration outside of the energy range returns zero."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    # Test below the energy range
    min_energy = energy[0] - 2
    max_energy = energy[0] - 1
    assert np.isclose(spec.integrate(min_energy, max_energy), 0.0), (
        "Integration below energy range should be zero"
    )

    # Test above the energy range
    min_energy = energy[-1] + 1
    max_energy = energy[-1] + 2
    assert np.isclose(spec.integrate(min_energy, max_energy), 0.0), (
        "Integration above energy range should be zero"
    )


def test_integrate_invalid_range() -> None:
    """Test that invalid integration ranges raise ValueError."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    with pytest.raises(
        ValueError,
        match="upper_energy must be greater than lower_energy",
    ):
        _ = spec.integrate(2, 2)

    with pytest.raises(
        ValueError,
        match="upper_energy must be greater than lower_energy",
    ):
        _ = spec.integrate(3, 1)


def test_write_csv() -> None:
    """Test writing and reading back spectrum CSV files."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    # Write the spectrum to a CSV
    filename = Path("test_spectrum.csv")
    if filename.exists():
        filename.unlink()
    spec.write_csv(filename)

    # Check that the file was created
    assert filename.exists(), "CSV file was not created"

    # Check that the contents are correct
    data_loaded = np.loadtxt(
        filename,
        delimiter=",",
        skiprows=2,
    )
    energy_lower_loaded = data_loaded[:, 0]
    energy_upper_loaded = data_loaded[:, 1]
    energy_loaded = np.concatenate([energy_lower_loaded, [energy_upper_loaded[-1]]])
    flux_loaded = data_loaded[:, 2]
    errors_loaded = data_loaded[:, 3]
    assert np.array_equal(energy, energy_loaded), (
        "Loaded energy edges do not match original"
    )
    assert np.array_equal(flux[:-1], flux_loaded), "Loaded flux does not match original"
    assert np.array_equal(errors[:-1], errors_loaded), (
        "Loaded errors do not match original"
    )

    # Clean up the file after testing
    filename.unlink()


def test_write_csv_filenames() -> None:
    """Test that write_csv defaults to correct file names."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1])

    # Test that file extension is added if not provided
    filename = Path("test_spectrum")
    expected_filename = filename.with_suffix(".csv")
    if expected_filename.exists():
        expected_filename.unlink()
    spec.write_csv(filename)
    assert expected_filename.exists(), "CSV file was not created"
    expected_filename.unlink()

    # Test with string filename not Path object
    filename = "test_spectrum"
    expected_filename = Path(filename + ".csv")
    if expected_filename.exists():
        expected_filename.unlink()
    spec.write_csv(filename)
    assert expected_filename.exists(), "CSV file was not created with string input"
    expected_filename.unlink()

    # Test that write_csv uses default filename derived from spectrum name.
    expected_filename = Path(f"{spec.name}.csv")
    if expected_filename.exists():
        expected_filename.unlink()
    spec.write_csv()  # no filename provided, should use default
    assert expected_filename.exists(), "write_csv() should create default filename"
    expected_filename.unlink()


def test_from_file() -> None:
    """Test creating a Spectrum from a CSV file."""
    energy, flux, errors = _mock_data()
    spec = Spectrum(energy=energy, flux=flux[:-1], errors=errors[:-1], name="test")

    # Write the spectrum to a CSV
    filename = Path("test_spectrum.csv")
    if filename.exists():
        filename.unlink()
    spec.write_csv(filename)

    # Load the spectrum back from the file
    loaded_spec = Spectrum.from_file(filename)

    # Test that the loaded spectrum matches the original
    assert isinstance(loaded_spec, Spectrum), "Loaded object is not a Spectrum"
    assert loaded_spec.name == spec.name, "Loaded spectrum name does not match"
    assert np.array_equal(loaded_spec.energy, spec.energy), (
        "Loaded energy edges do not match original"
    )
    assert np.array_equal(loaded_spec.flux, spec.flux), (
        "Loaded flux does not match original"
    )
    assert np.array_equal(loaded_spec.errors, spec.errors), (
        "Loaded errors do not match original"
    )

    # Clean up the file after testing
    filename.unlink()
