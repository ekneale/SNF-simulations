"""Unit tests for loading antineutrino spectra data."""

import numpy as np
import ROOT
from snf_simulations.load_and_scale import load_equal_scaled
from snf_simulations.load_data import load_antineutrino_data
from snf_simulations.load_spec import load_equal, load_spec
from snf_simulations.scale import scale


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


def _test_load_spec(isotope_data, isotope_name="test_isotope"):
    """Test that spectra histograms can be created and have the correct properties."""
    # NOTE: This will produce ROOT windows when run, which may need to be closed manually.
    #       Run pytest with -s to allow closing each window as they pop up.
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
    # Since we're using the energy values as bin edges, there are len(energy)-1 bins.
    expected_lowers = energy[:-1]
    expected_uppers = energy[1:]
    expected_contents = dN[:-1]
    expected_errors = errors[:-1]
    for i, nbin in enumerate(range(1, spec.GetNbinsX() + 1)):
        # Check that bin edges, content and errors match the expected values
        assert spec.GetBinLowEdge(nbin) == expected_lowers[i], (
            f"{isotope_name} spectrum bin {nbin} lower edge mismatch"
        )
        assert (
            spec.GetBinLowEdge(nbin) + spec.GetBinWidth(nbin) == expected_uppers[i]
        ), f"{isotope_name} spectrum bin {nbin} upper edge mismatch"
        assert spec.GetBinContent(nbin) == expected_contents[i], (
            f"{isotope_name} spectrum bin {nbin} content mismatch"
        )
        assert spec.GetBinError(nbin) == expected_errors[i], (
            f"{isotope_name} spectrum bin {nbin} error mismatch"
        )

    return spec


def test_load_spec_mock():
    """Test that spectra histograms can be created with mock data."""
    # Create some fake data
    # NOTE: Need to set dtype to float to keep ROOT happy
    energy = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3], dtype=float)
    dN = np.array([10, 20, 30, 40, 50, 60, 70], dtype=float)
    errors = np.array([1, 2, 3, 4, 5, 6, 7], dtype=float)
    data = np.column_stack((energy, dN, errors))

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
    # (since the bin widths don't have to be equal - see load_equal).
    spec = _test_load_spec(data)
    centres = np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75])
    content = np.array([10, 20, 30, 40, 50, 60])
    errors = np.array([1, 2, 3, 4, 5, 6])
    for i, nbin in enumerate(range(1, spec.GetNbinsX() + 1)):
        assert spec.GetBinCenter(nbin) == centres[i], (
            f"test_isotope spectrum bin {nbin} centre mismatch"
        )
        assert spec.GetBinContent(nbin) == content[i], (
            f"test_isotope spectrum bin {nbin} content mismatch"
        )
        assert spec.GetBinError(nbin) == errors[i], (
            f"test_isotope spectrum bin {nbin} error mismatch"
        )


def test_load_spec_real():
    """Test that spectra histograms can be created for the included isotope data."""
    data = load_antineutrino_data()
    for i, isotope_data in enumerate(data):
        isotope_name = f"isotope_{i}"
        _test_load_spec(isotope_data, isotope_name)


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


def _test_load_equal(isotope_data, isotope_name="test_isotope", min_E=0, max_E=None):
    """Test that spectra can be loaded with equal bin widths."""
    energy = isotope_data[:, 0]
    dN = isotope_data[:, 1]
    errors = isotope_data[:, 2]
    if max_E is None:
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
    expected_edges = np.linspace(min_E, max_E, (max_E - min_E) + 1)
    expected_centres = expected_edges[:-1] + 0.5
    for i, nbin in enumerate(range(1, spec.GetNbinsX() + 1)):
        # Check that bin widths are equally spaced integers
        assert spec.GetBinWidth(nbin) == 1, (
            f"{isotope_name} spectrum bin {nbin} width is not 1"
        )
        # Check that bin centers match expected values
        assert spec.GetBinCenter(nbin) == expected_centres[i], (
            f"{isotope_name} spectrum bin {nbin} center mismatch"
        )
        # Check that bin edges match expected values
        assert spec.GetBinLowEdge(nbin) == expected_edges[i], (
            f"{isotope_name} spectrum bin {nbin} lower edge mismatch"
        )

    # Test bin contents and errors
    original_centres = energy[:-1] + np.diff(energy) / 2  # energy values are bin edges
    original_content = dN[:-1]
    original_errors = errors[:-1]
    expected_contents, expected_errors = _linear_interpolate_with_errors(
        original_centres, original_content, original_errors, expected_centres
    )
    for i, nbin in enumerate(range(1, spec.GetNbinsX() + 1)):
        bin_center = spec.GetBinCenter(nbin)
        bin_content = spec.GetBinContent(nbin)
        bin_error = spec.GetBinError(nbin)

        assert np.isclose(bin_center, expected_centres[i]), (
            f"{isotope_name} spectrum bin {nbin} center mismatch"
        )
        assert np.isclose(bin_content, expected_contents[i]), (
            f"{isotope_name} spectrum bin {nbin} content mismatch"
        )
        assert np.isclose(bin_error, expected_errors[i]), (
            f"{isotope_name} spectrum bin {nbin} error mismatch"
        )

    return spec


def test_load_equal_mock():
    """Test that equal spectra histograms can be created with mock data."""
    # Create some fake data
    energy = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3], dtype=float)
    dN = np.array([10, 20, 30, 40, 50, 60, 70], dtype=float)
    errors = np.array([1, 2, 3, 4, 5, 6, 7], dtype=float)
    data = np.column_stack((energy, dN, errors))

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
    # The contents are found at the new centres by linear interpolation, while the errors are
    # calculated in quadrature based on the closest bins.
    # e.g. for new bin 1 (centre 0.5):
    #  lower old bin 1 centre = 0.25, error = 1.0
    #  upper old bin 2 centre = 0.75, error = 2.0
    #  distance to lower = 0.25
    #  distance to upper = 0.25
    #  total distance = 0.5
    #  weight lower = 0.25 / 0.5 = 0.5
    #  weight upper = 0.25 / 0.5 = 0.5
    #  new error = sqrt(0.5^2 * 1.0^2 + 0.5^2 * 2.0^2) = 1.118...
    spec = _test_load_equal(data)
    centres = np.array([0.5, 1.5, 2.5])
    content = np.array([15.0, 35.0, 55.0])
    errors = np.array(
        [
            np.sqrt((0.25 / 0.5) ** 2 * 1.0**2 + (0.25 / 0.5) ** 2 * 2.0**2),
            np.sqrt((0.25 / 0.5) ** 2 * 3.0**2 + (0.25 / 0.5) ** 2 * 4.0**2),
            np.sqrt((0.25 / 0.5) ** 2 * 5.0**2 + (0.25 / 0.5) ** 2 * 6.0**2),
        ]
    )
    for i, nbin in enumerate(range(1, spec.GetNbinsX() + 1)):
        assert spec.GetBinCenter(nbin) == centres[i], (
            f"test_isotope spectrum bin {nbin} centre mismatch"
        )
        assert spec.GetBinContent(nbin) == content[i], (
            f"test_isotope spectrum bin {nbin} content mismatch"
        )
        assert np.isclose(spec.GetBinError(nbin), errors[i]), (
            f"test_isotope spectrum bin {nbin} error mismatch"
        )


def test_load_equal_mock_extrapolate():
    """Test that equal spectra histograms can be extrapolated beyond the original range."""
    # Create some fake data
    energy = np.array([0, 2, 4, 6], dtype=float)
    dN = np.array([10, 20, 30, 40], dtype=float)
    errors = np.array([1, 2, 3, 4], dtype=float)
    data = np.column_stack((energy, dN, errors))

    # The load_spec function will create bins at every integer from min_E to max_E.
    min_E = -1
    max_E = 7
    # The test output should go from this (from load_spec):
    # bin  lower  centre  upper  content  error
    #   1    0.0     1.0    2.0     10.0    1.0
    #   2    2.0     3.0    4.0     20.0    2.0
    #   3    4.0     5.0    6.0     30.0    3.0
    # to this:
    # bin  lower  centre  upper  content  error
    #   1   -1.0    -0.5    0.0     10.0    1.0
    #   2    0.0     0.5    1.0     10.0    1.0
    #   3    1.0     1.5    2.0     12.5    0.9013...
    #   4    2.0     2.5    3.0     17.5    1.5206...
    #   5    3.0     3.5    4.0     22.5    1.6770...
    #   6    4.0     4.5    5.0     27.5    2.3048...
    #   7    5.0     5.5    6.0     30.0    3.0
    #   8    6.0     6.5    7.0     30.0    3.0
    # NOTE: Any bins below the lowest original bin centre (1.0) or above the highest (5.0)
    #       just take the content and error of the first/last original bin.
    spec = _test_load_equal(data, min_E=min_E, max_E=max_E)
    centres = np.array([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
    content = np.array([10.0, 10.0, 12.5, 17.5, 22.5, 27.5, 30.0, 30.0])
    errors = np.array(
        [
            1.0,
            1.0,
            np.sqrt((1.5 / 2.0) ** 2 * 1.0**2 + (0.5 / 2.0) ** 2 * 2.0**2),
            np.sqrt((0.5 / 2.0) ** 2 * 1.0**2 + (1.5 / 2.0) ** 2 * 2.0**2),
            np.sqrt((1.5 / 2.0) ** 2 * 2.0**2 + (0.5 / 2.0) ** 2 * 3.0**2),
            np.sqrt((0.5 / 2.0) ** 2 * 2.0**2 + (1.5 / 2.0) ** 2 * 3.0**2),
            3.0,
            3.0,
        ]
    )
    for i, nbin in enumerate(range(1, spec.GetNbinsX() + 1)):
        assert spec.GetBinCenter(nbin) == centres[i], (
            f"test_isotope spectrum bin {nbin} centre mismatch"
        )
        assert spec.GetBinContent(nbin) == content[i], (
            f"test_isotope spectrum bin {nbin} content mismatch"
        )
        assert np.isclose(spec.GetBinError(nbin), errors[i]), (
            f"test_isotope spectrum bin {nbin} error mismatch"
        )


def test_load_equal_real():
    """Test that equal spectra histograms can be created for the included isotope data."""
    data = load_antineutrino_data()
    for i, isotope_data in enumerate(data):
        isotope_name = f"isotope_{i}"
        _test_load_equal(isotope_data, isotope_name=isotope_name)


def _isotope_activity(mass_kg, atomic_mass, half_life_yrs, removal_time_yrs):
    """Calculate the activity of an isotope spectrum after a given time."""
    # NOTE: Taken from scale.scale (and cleaned up a bit)
    # TODO: Move this to scale.py and reuse in both places,
    #       Then add unit tests for this function too.
    # Convert mass to number of atoms (kg to g, then to moles, then to atoms)
    N0 = (mass_kg * 1000 / atomic_mass) * 6.022e23
    # Calculate initial activity (in decays per second aka Becquerels)
    lambda_ = np.log(2) / (half_life_yrs * 365 * 24 * 60 * 60)
    A0 = N0 * lambda_
    # Calculate activity after removal time (still in decays per second)
    A = A0 * np.exp(-1 * lambda_ * removal_time_yrs * 365 * 24 * 60 * 60)
    return A


def _test_scale(spec, m, mr, half_life_yrs, removal_time_yrs):
    """Test that scaling spectra works as expected."""
    # Scale the spectrum
    # NOTE ROOT edits in place, so make a copy to keep the original unscaled
    scaled_spec = scale(spec.Clone(), m, mr, half_life_yrs, removal_time_yrs)

    # Test basic properties
    assert isinstance(scaled_spec, ROOT.TH1D), "Scaled spectrum is not a TH1D"
    assert scaled_spec.GetNbinsX() == spec.GetNbinsX(), (
        "Scaled spectrum has different number of bins"
    )

    # Test that each bin content and error has been scaled correctly
    expected_activity = _isotope_activity(m, mr, half_life_yrs, removal_time_yrs)
    for nbin in range(1, spec.GetNbinsX() + 1):
        original_content = spec.GetBinContent(nbin)
        scaled_content = scaled_spec.GetBinContent(nbin)
        expected_content = original_content * expected_activity
        assert np.isclose(scaled_content, expected_content), (
            f"Scaled spectrum bin {nbin} content mismatch"
        )

        original_error = spec.GetBinError(nbin)
        scaled_error = scaled_spec.GetBinError(nbin)
        expected_error = original_error * expected_activity
        assert np.isclose(scaled_error, expected_error), (
            f"Scaled spectrum bin {nbin} error mismatch"
        )

    return scaled_spec


def test_scale_mock():
    """Test that scaling works on a mock spectrum."""
    # Create a fake spectrum
    energy = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3], dtype=float)
    dN = np.array([10, 20, 30, 40, 50, 60, 70], dtype=float)
    errors = np.array([1, 2, 3, 4, 5, 6, 7], dtype=float)
    spec = load_spec(
        Energy=energy,
        dN=dN,
        errors=errors,
        isotope='test_isotope',
    )

    # Define scaling parameters
    m = 1.0  # kg
    mr = 100.0  # g/mol
    half_life_yrs = 1.0  # year
    removal_time_yrs = 0.5  # years

    scaled_spec = _test_scale(spec, m, mr, half_life_yrs, removal_time_yrs)
    A = 9.359326705935426e16  # Pre-calculated expected activity for these parameters
    expected_contents = dN[:-1] * A
    expected_errors = errors[:-1] * A
    for i, nbin in enumerate(range(1, scaled_spec.GetNbinsX() + 1)):
        assert np.isclose(scaled_spec.GetBinContent(nbin), expected_contents[i]), (
            f"Scaled mock spectrum bin {nbin} content mismatch"
        )
        assert np.isclose(scaled_spec.GetBinError(nbin), expected_errors[i]), (
            f"Scaled mock spectrum bin {nbin} error mismatch"
        )


def _test_load_and_scale(
    data, name, m, mr, half_life_yrs, removal_time_yrs, max_E, min_E=0
):
    """Test that loading and scaling spectra works as expected."""
    # Use the combined function
    ls_spec = load_equal_scaled(
        data, max_E, name, name, m, mr, half_life_yrs, removal_time_yrs, min_E=min_E
    )
    # Use the individual functions for comparison
    loaded_spec = load_equal(
        name=name,
        isotope=name,
        E=data[:, 0],
        dN=data[:, 1],
        error=data[:, 2],
        max_E=max_E,
        min_E=min_E,
    )
    scaled_spec = scale(loaded_spec.Clone(), m, mr, half_life_yrs, removal_time_yrs)

    # Test that both methods give the same result
    assert ls_spec.GetNbinsX() == scaled_spec.GetNbinsX(), (
        "Load and scale spectrum has different number of bins"
    )
    for nbin in range(1, ls_spec.GetNbinsX() + 1):
        ls_content = ls_spec.GetBinContent(nbin)
        scaled_content = scaled_spec.GetBinContent(nbin)
        assert np.isclose(ls_content, scaled_content), (
            f"Load and scale spectrum bin {nbin} content mismatch"
        )

        ls_error = ls_spec.GetBinError(nbin)
        scaled_error = scaled_spec.GetBinError(nbin)
        assert np.isclose(ls_error, scaled_error), (
            f"Load and scale spectrum bin {nbin} error mismatch"
        )


def test_load_and_scale_mock():
    """Test that loading and scaling works on a mock spectrum."""
    # Create a fake spectrum
    energy = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3], dtype=float)
    dN = np.array([10, 20, 30, 40, 50, 60, 70], dtype=float)
    errors = np.array([1, 2, 3, 4, 5, 6, 7], dtype=float)
    data = np.column_stack((energy, dN, errors))

    # Define scaling parameters
    m = 1.0  # kg
    mr = 100.0  # g/mol
    half_life_yrs = 1.0  # year
    removal_time_yrs = 0.5  # years

    _test_load_and_scale(
        data,
        "test_isotope",
        m,
        mr,
        half_life_yrs,
        removal_time_yrs,
        max_E=3,
    )
