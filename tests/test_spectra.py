"""Unit tests for loading antineutrino spectra data."""

import numpy as np
import ROOT
from snf_simulations.load_data import load_antineutrino_data
from snf_simulations.load_spec import load_equal, load_spec


def test_load_antineutrino_data():
    """Test that antineutrino data is imported correctly."""
    data = load_antineutrino_data()
    assert isinstance(data, tuple), "Loaded data is not a tuple"
    assert len(data) == 16, f"Loaded {len(data)} isotopes, expected 16"
    for i, isotope_data in enumerate(data):
        # Basic structure checks
        assert isinstance(isotope_data, np.ndarray), (
            f"Isotope {i} data is not a numpy array"
        )
        assert isotope_data.shape[1] == 3, (
            f"Isotope {i} data has {isotope_data.shape[1]} columns, expected 3"
        )
        assert isotope_data.shape[0] > 0, f"Isotope {i} data has no rows"

        # Physical value checks
        assert np.all(isotope_data[:, 0] >= 0), f"Isotope {i} has negative energies"
        assert np.all(isotope_data[:, 1] >= 0), f"Isotope {i} has negative flux values"
        assert np.all(isotope_data[:, 2] >= 0), f"Isotope {i} has negative error values"


def test_load_spec(isotope_data=None, isotope_name="test_isotope"):
    """Test that spectra histograms can be created and have the correct properties."""
    # NOTE: This will produce ROOT windows when run, which may need to be closed manually.
    #       Run pytest with -s to allow closing each window as they pop up.
    if isotope_data is None:
        # Create some fake data
        # NOTE: Need to set dtype to float to keep ROOT happy
        energy = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3], dtype=float)
        dN = np.array([10, 20, 30, 40, 50, 60, 70], dtype=float)
        errors = np.array([1, 2, 3, 4, 5, 6, 7], dtype=float)
    else:
        energy = isotope_data[:, 0]
        dN = isotope_data[:, 1]
        errors = isotope_data[:, 2]

    # Create the spectrum
    spec = load_spec(
        Energy=energy,
        dN=dN,
        errors=errors,
        isotope=isotope_name,
    )

    # Test basic properties
    assert isinstance(spec, ROOT.TH1D), f"{isotope_name} spectrum is not a TH1D"
    assert spec.GetName() == isotope_name, f"{isotope_name} spectrum name mismatch"
    assert spec.GetNbinsX() == len(energy) - 1, (
        f"{isotope_name} spectrum has {spec.GetNbinsX()} bins instead of {len(energy) - 1}"
    )
    assert spec.GetNbinsY() == 1, f"{isotope_name} spectrum has more than 1 Y bin"
    assert spec.GetEntries() == len(energy), f"{isotope_name} spectrum entries mismatch"

    # Test bin contents and errors
    # ROOT bins are 1-indexed, with bin 0 being underflow and bin N+1 overflow.
    # Each bin is defined with the input energy values being the edges,
    # so the lower edge of bin 1 is energy[0], the upper edge of bin 1 (and the lower edge
    # of bin 2) is energy[1], etc.
    # The test output should look like this:
    # bin  lower  centre  upper  content  error
    #   1    0.0    0.25    0.5     10.0    1.0
    #   2    0.5    0.75    1.0     20.0    2.0
    #   3    1.0    1.25    1.5     30.0    3.0
    #   4    1.5    1.75    2.0     40.0    4.0
    #   5    2.0    2.25    2.5     50.0    5.0
    #   6    2.5    2.75    3.0     60.0    6.0
    # NOTE: We lose the last data point (3, 70, 7), as it defines the upper edge of the final bin.
    # If we had a bin 7 then we wouldn't know what the upper limit would be
    # (since the bin widths don't have to be equal - see load_equal)
    expected_lowers = energy[:-1]
    expected_uppers = energy[1:]
    expected_contents = dN[:-1]
    expected_errors = errors[:-1]
    for i in range(1, spec.GetNbinsX() + 1):
        # Check that bin edges, content and errors match the expected values
        assert spec.GetBinLowEdge(i) == expected_lowers[i - 1], (
            f"{isotope_name} spectrum bin {i} lower edge mismatch"
        )
        assert spec.GetBinLowEdge(i) + spec.GetBinWidth(i) == expected_uppers[i - 1], (
            f"{isotope_name} spectrum bin {i} upper edge mismatch"
        )
        assert spec.GetBinContent(i) == expected_contents[i - 1], (
            f"{isotope_name} spectrum bin {i} content mismatch"
        )
        assert spec.GetBinError(i) == expected_errors[i - 1], (
            f"{isotope_name} spectrum bin {i} error mismatch"
        )


def test_load_spec_with_data():
    """Test that spectra histograms can be created for the included isotope data."""
    data = load_antineutrino_data()
    for i, isotope_data in enumerate(data):
        isotope_name = f"isotope_{i}"
        test_load_spec(isotope_data=isotope_data, isotope_name=isotope_name)


def _linear_interpolate_with_errors(
    original_centres, original_content, original_errors, new_centres
):
    """Linearly interpolate content and propagate errors to new bin centers."""
    # TODO: Eventually this should replace the ROOT code in load_spec.py

    # Interpolate new content values
    new_content = np.interp(new_centres, original_centres, original_content)

    # Propagate errors
    new_errors = np.zeros_like(new_centres)
    for i, centre in enumerate(new_centres):
        # Check for extrapolation cases
        if centre < original_centres[0]:
            # Extrapolation below the first bin centre - use first error
            new_errors[i] = original_errors[0]
            continue
        elif centre >= original_centres[-1]:
            # Extrapolation above the last bin centre - use last error
            new_errors[i] = original_errors[-1]
            continue

        # Find the two closest original bin centers.
        # Using np.searchsorted finds the "insertion point" for the new centre,
        # i.e. the index of where it would go to keep the array sorted.
        # So if the original_centres are [1, 2, 3] and centre is 2.5, idx will be 2
        # as it would fit between 2 (index 1) and 3 (index 2).
        # Therefore the surrounding bins are at idx-1 and idx.
        idx = np.searchsorted(original_centres, centre)
        lower_idx = idx - 1
        upper_idx = idx

        # Calculate new error by propagating errors from the two surrounding bins,
        # weighted by distance to the new centre.
        c_lower = original_centres[lower_idx]
        c_upper = original_centres[upper_idx]
        err_lower = original_errors[lower_idx]
        err_upper = original_errors[upper_idx]
        weight_lower = (c_upper - centre) / (c_upper - c_lower)
        weight_upper = (centre - c_lower) / (c_upper - c_lower)
        new_errors[i] = np.sqrt(
            weight_lower**2 * err_lower**2 + weight_upper**2 * err_upper**2
        )

    return new_content, new_errors


def test_load_equal(isotope_data=None, isotope_name="test_isotope"):
    """Test that spectra can be loaded with equal bin widths."""
    if isotope_data is None:
        # Create some fake data
        # NOTE: Need to set dtype to float to keep ROOT happy
        energy = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3], dtype=float)
        dN = np.array([10, 20, 30, 40, 50, 60, 70], dtype=float)
        errors = np.array([1, 2, 3, 4, 5, 6, 7], dtype=float)
    else:
        energy = isotope_data[:, 0]
        dN = isotope_data[:, 1]
        errors = isotope_data[:, 2]
    min_E = 0  # BUG: The function forces the lowest bin to 0.5, so min_E has to be 0.
    max_E = int(np.max(energy))  # rounds down to nearest int

    # Create the spectrum
    spec = load_equal(
        name=isotope_name,
        isotope=isotope_name,
        E=energy,
        dN=dN,
        error=errors,
        max_E=max_E,
        min_E=min_E,
    )

    # Test basic properties
    assert isinstance(spec, ROOT.TH1D), f"{isotope_name} spectrum is not a TH1D"
    assert spec.GetName() == isotope_name, f"{isotope_name} spectrum name mismatch"
    assert spec.GetNbinsX() == max_E - min_E, (
        f"{isotope_name} spectrum has {spec.GetNbinsX()} bins instead of {max_E - min_E}"
    )
    assert spec.GetNbinsY() == 1, f"{isotope_name} spectrum has more than 1 Y bin"
    assert spec.GetEntries() == max_E - min_E, (
        f"{isotope_name} spectrum entries mismatch"
    )

    # Test bin spacing
    # The bins should now be equally spaced integers from min_E to max_E.
    for i in range(1, max_E + 1):
        # Check that bin widths are equally spaced integers
        assert spec.GetBinLowEdge(i) == (i - 1), (
            f"{isotope_name} spectrum bin {i} lower edge is not {i - 1}"
        )
        assert spec.GetBinWidth(i) == 1, (
            f"{isotope_name} spectrum bin {i} width is not 1"
        )

    # Test bin contents and errors
    # The test output should go from this (from load_spec):
    # bin  lower  centre  upper  content  error
    #   1    0.0    0.25    0.5     10.0    1.0
    #   2    0.5    0.75    1.0     20.0    2.0
    #   3    1.0    1.25    1.5     30.0    3.0
    #   4    1.5    1.75    2.0     40.0    4.0
    #   5    2.0    2.25    2.5     50.0    5.0
    #   6    2.5    2.75    3.0     60.0    6.0
    # to this:
    # bin  lower  centre  upper  content  error
    #   1    0.0     0.5    1.0     15.0    1.118...
    #   2    1.0     1.5    2.0     35.0    2.5
    #   3    2.0     2.5    3.0     55.0    3.905...
    # The contents are found at the new centres by linear interpolation,
    # using the ROOT TH1D::Interpolate method.
    # The errors are a bit more complicated, they're calculated in quadrature based on the
    # closest bins.
    # e.g. for new bin 1 (centre 0.5):
    #  lower old bin 1 centre = 0.25, error = 1.0
    #  upper old bin 2 centre = 0.75, error = 2.0
    #  distance to lower = 0.25
    #  distance to upper = 0.25
    #  total distance = 0.5
    #  weight lower = 0.25 / 0.5 = 0.5
    #  weight upper = 0.25 / 0.5 = 0.5
    #  new error = sqrt(0.5^2 * 1.0^2 + 0.5^2 * 2.0^2) = 1.118...
    original_centres = energy[:-1] + np.diff(energy) / 2
    original_content = dN[:-1]
    original_errors = errors[:-1]
    expected_centres = np.linspace(0.5, max_E - 0.5, max_E)
    expected_contents, expected_errors = _linear_interpolate_with_errors(
        original_centres, original_content, original_errors, expected_centres
    )
    for i in range(1, max_E + 1):
        bin_center = spec.GetBinCenter(i)
        bin_content = spec.GetBinContent(i)
        bin_error = spec.GetBinError(i)

        assert np.isclose(bin_center, expected_centres[i - 1]), (
            f"{isotope_name} spectrum test bin {i} center mismatch"
        )
        assert np.isclose(bin_content, expected_contents[i - 1]), (
            f"{isotope_name} spectrum test bin {i} content mismatch"
        )
        assert np.isclose(bin_error, expected_errors[i - 1]), (
            f"{isotope_name} spectrum test bin {i} error mismatch"
        )


def test_load_equal_with_data():
    """Test that equal spectra histograms can be created for the included isotope data."""
    data = load_antineutrino_data()
    for i, isotope_data in enumerate(data):
        isotope_name = f"isotope_{i}"
        test_load_equal(isotope_data=isotope_data, isotope_name=isotope_name)
