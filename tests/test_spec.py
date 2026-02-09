"""Unit tests for loading antineutrino spectra data."""

from pathlib import Path

import numpy as np
import ROOT

from snf_simulations.data import load_antineutrino_data, load_isotope_data
from snf_simulations.physics import get_isotope_activity
from snf_simulations.spec import (
    add_spec,
    create_spec,
    equalise_spec,
    load_spec,
    sample_spec,
    scale_spec,
    write_spec,
)

from .test_commandline import _spec_to_arrays

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers


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


def _test_create_spec(
    isotope_data: np.ndarray,
    isotope_name: str = "test_isotope",
) -> ROOT.TH1D:
    """Test that spectra histograms can be created and have the correct properties."""
    # Create the spectrum
    energy = isotope_data[:, 0]
    dn_de = isotope_data[:, 1]
    errors = isotope_data[:, 2]
    spec = create_spec(
        energy=energy,
        dn_de=dn_de,
        errors=errors,
        name=isotope_name,
    )

    # Test basic properties
    assert isinstance(spec, ROOT.TH1D), f"{isotope_name} spectrum is not a TH1D"
    assert spec.GetName() == isotope_name, f"{isotope_name} spectrum name mismatch"
    assert spec.GetNbinsX() == len(energy) - 1, (
        f"{isotope_name} spectrum has {spec.GetNbinsX()} bins"
        f" instead of {len(energy) - 1}"
    )
    assert spec.GetNbinsY() == 1, f"{isotope_name} spectrum has more than 1 Y bin"
    assert spec.GetEntries() == len(energy), f"{isotope_name} spectrum entries mismatch"

    # Test bin contents and errors
    # ROOT bins are 1-indexed, with bin 0 being underflow and bin N+1 overflow.
    # Each bin is defined with the input energy values being the edges,
    # so the lower edge of bin 1 is energy[0], the upper edge of bin 1 (and the
    # lower edge of bin 2) is energy[1], etc.
    # Since we're using the energy values as bin edges, there are len(energy)-1 bins.
    expected_lowers = energy[:-1]
    expected_uppers = energy[1:]
    expected_contents = dn_de[:-1]
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


def test_create_spec_mock() -> None:
    """Test that spectra histograms can be created with mock data."""
    # Create some fake data
    # NOTE: Need to set dtype to float to keep ROOT happy
    energy = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3], dtype=float)
    dn_de = np.array([10, 20, 30, 40, 50, 60, 70], dtype=float)
    errors = np.array([1, 2, 3, 4, 5, 6, 7], dtype=float)
    data = np.column_stack((energy, dn_de, errors))

    # The test output should look like this:
    # bin  lower  centre  upper  content  error
    #   1    0.0    0.25    0.5     10.0    1.0
    #   2    0.5    0.75    1.0     20.0    2.0
    #   3    1.0    1.25    1.5     30.0    3.0
    #   4    1.5    1.75    2.0     40.0    4.0
    #   5    2.0    2.25    2.5     50.0    5.0
    #   6    2.5    2.75    3.0     60.0    6.0
    # NOTE: We lose the last data point (3, 70, 7), as it defines the upper edge of
    # the final bin.
    # If we had a bin 7 then we wouldn't know what the upper limit would be
    # (since the bin widths don't have to be equal - see equalise_spec).
    spec = _test_create_spec(data)
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


def test_create_spec_real() -> None:
    """Test that spectra histograms can be created for the included isotope data."""
    molar_masses, _ = load_isotope_data()
    isotopes = list(molar_masses.keys())
    data = load_antineutrino_data(isotopes)
    for isotope, isotope_data in data.items():
        _test_create_spec(isotope_data, isotope_name=isotope)


def _linear_interpolate_with_errors(
    original_centres: np.ndarray,
    original_content: np.ndarray,
    original_errors: np.ndarray,
    new_centres: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Linearly interpolate content and propagate errors to new bin centers."""
    # TODO: Eventually this should replace the ROOT code in spec.equalise_spec

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
        if centre >= original_centres[-1]:
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
            weight_lower**2 * err_lower**2 + weight_upper**2 * err_upper**2,
        )

    return new_content, new_errors


def _test_equalise_spec(
    isotope_data: np.ndarray,
    isotope_name: str = "test_isotope",
    min_energy: int = 0,
    max_energy: int | None = None,
) -> ROOT.TH1D:
    """Test that spectra can be loaded with equal bin widths."""
    energy = isotope_data[:, 0]
    dn_de = isotope_data[:, 1]
    errors = isotope_data[:, 2]
    if max_energy is None:
        max_energy = int(np.max(energy))  # rounds down to nearest int

    # Create the spectrum and equalise it
    spec = create_spec(
        energy=energy,
        dn_de=dn_de,
        errors=errors,
        name=isotope_name,
    )
    spec = equalise_spec(spec, max_energy, min_energy)

    # Test basic properties
    assert isinstance(spec, ROOT.TH1D), f"{isotope_name} spectrum is not a TH1D"
    assert spec.GetName() == isotope_name, f"{isotope_name} spectrum name mismatch"
    assert spec.GetNbinsX() == max_energy - min_energy, (
        f"{isotope_name} spectrum has {spec.GetNbinsX()} bins"
        f" instead of {max_energy - min_energy}"
    )
    assert spec.GetNbinsY() == 1, f"{isotope_name} spectrum has more than 1 Y bin"
    assert spec.GetEntries() == max_energy - min_energy, (
        f"{isotope_name} spectrum entries mismatch"
    )

    # Test bin spacing
    # The bins should now be equally spaced integers from min_energy to max_energy.
    expected_edges = np.linspace(min_energy, max_energy, (max_energy - min_energy) + 1)
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
    original_content = dn_de[:-1]
    original_errors = errors[:-1]
    expected_contents, expected_errors = _linear_interpolate_with_errors(
        original_centres,
        original_content,
        original_errors,
        expected_centres,
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


def test_equalise_spec_mock() -> None:
    """Test that equal spectra histograms can be created with mock data."""
    # Create some fake data
    energy = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3], dtype=float)
    dn_de = np.array([10, 20, 30, 40, 50, 60, 70], dtype=float)
    errors = np.array([1, 2, 3, 4, 5, 6, 7], dtype=float)
    data = np.column_stack((energy, dn_de, errors))

    # The test output should go from this (from create_spec):
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
    # while the errors are calculated in quadrature based on the closest bins.
    # e.g. for new bin 1 (centre 0.5):
    #  lower old bin 1 centre = 0.25, error = 1.0
    #  upper old bin 2 centre = 0.75, error = 2.0
    #  distance to lower = 0.25
    #  distance to upper = 0.25
    #  total distance = 0.5
    #  weight lower = 0.25 / 0.5 = 0.5
    #  weight upper = 0.25 / 0.5 = 0.5
    #  new error = sqrt(0.5^2 * 1.0^2 + 0.5^2 * 2.0^2) = 1.118...
    spec = _test_equalise_spec(data)
    centres = np.array([0.5, 1.5, 2.5])
    content = np.array([15.0, 35.0, 55.0])
    errors = np.array(
        [
            np.sqrt((0.25 / 0.5) ** 2 * 1.0**2 + (0.25 / 0.5) ** 2 * 2.0**2),
            np.sqrt((0.25 / 0.5) ** 2 * 3.0**2 + (0.25 / 0.5) ** 2 * 4.0**2),
            np.sqrt((0.25 / 0.5) ** 2 * 5.0**2 + (0.25 / 0.5) ** 2 * 6.0**2),
        ],
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


def test_equalise_spec_mock_extrapolate() -> None:
    """Test that spectra histograms can be extrapolated beyond the original range."""
    # Create some fake data
    energy = np.array([0, 2, 4, 6], dtype=float)
    dn_de = np.array([10, 20, 30, 40], dtype=float)
    errors = np.array([1, 2, 3, 4], dtype=float)
    data = np.column_stack((energy, dn_de, errors))

    # The create_spec function will create bins at every integer from min_energy
    # to max_energy.
    min_energy = -1
    max_energy = 7
    # The test output should go from this (from create_spec):
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
    # NOTE: Any bins below the lowest original bin centre (1.0)
    #   or above the highest (5.0) just take the content and error of the
    #   first/last original bin.
    spec = _test_equalise_spec(data, min_energy=min_energy, max_energy=max_energy)
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
        ],
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


def test_equalise_spec_real() -> None:
    """Test that equal histograms can be created for the included isotope data."""
    molar_masses, _ = load_isotope_data()
    isotopes = list(molar_masses.keys())
    data = load_antineutrino_data(isotopes)
    for isotope, isotope_data in data.items():
        _test_equalise_spec(isotope_data, isotope_name=isotope)


def _test_scale_spec(
    spec: ROOT.TH1D,
    mass: float,
    molar_mass: float,
    half_life: float,
    removal_time: float,
) -> ROOT.TH1D:
    """Test that scaling spectra works as expected."""
    # Scale the spectrum
    # NOTE ROOT edits in place, so make a copy to keep the original unscaled
    scaled_spec = scale_spec(spec.Clone(), mass, molar_mass, half_life, removal_time)

    # Test basic properties
    assert isinstance(scaled_spec, ROOT.TH1D), "Scaled spectrum is not a TH1D"
    assert scaled_spec.GetNbinsX() == spec.GetNbinsX(), (
        "Scaled spectrum has different number of bins"
    )

    # Test that each bin content and error has been scaled correctly
    expected_activity = get_isotope_activity(mass, molar_mass, half_life, removal_time)
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


def test_scale_mock() -> None:
    """Test that scaling works on a mock spectrum."""
    # Create a fake spectrum
    energy = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3], dtype=float)
    dn_de = np.array([10, 20, 30, 40, 50, 60, 70], dtype=float)
    errors = np.array([1, 2, 3, 4, 5, 6, 7], dtype=float)
    spec = create_spec(
        energy=energy,
        dn_de=dn_de,
        errors=errors,
        name="test_isotope",
    )

    # Define scaling parameters
    mass = 1.0  # kg
    molar_mass = 100.0  # g/mol
    half_life = 1.0  # year
    removal_time = 0.5  # years

    scaled_spec = _test_scale_spec(spec, mass, molar_mass, half_life, removal_time)
    activity = 9.359326705935426e16  # Pre-calculated activity for these parameters
    expected_contents = dn_de[:-1] * activity
    expected_errors = errors[:-1] * activity
    for i, nbin in enumerate(range(1, scaled_spec.GetNbinsX() + 1)):
        assert np.isclose(scaled_spec.GetBinContent(nbin), expected_contents[i]), (
            f"Scaled mock spectrum bin {nbin} content mismatch"
        )
        assert np.isclose(scaled_spec.GetBinError(nbin), expected_errors[i]), (
            f"Scaled mock spectrum bin {nbin} error mismatch"
        )


def _test_load_and_scale(  # noqa: PLR0913
    data: np.ndarray,
    name: str,
    mass: float,
    molar_mass: float,
    half_life: float,
    removal_time: float,
    max_energy: int,
    min_energy: int = 0,
) -> ROOT.TH1D:
    """Test that loading and scaling spectra works as expected."""
    # TODO: This has too many arguments but it mirrors load_spec, see it for details.
    # Use the combined function
    ls_spec = load_spec(
        data,
        name,
        mass,
        molar_mass,
        half_life,
        removal_time,
        max_energy=None,
        min_energy=0,
    )
    # Use the individual functions for comparison
    spec = create_spec(
        energy=data[:, 0],
        dn_de=data[:, 1],
        errors=data[:, 2],
        name=name,
    )
    equalised_spec = equalise_spec(spec, max_energy, min_energy)
    scaled_spec = scale_spec(
        equalised_spec.Clone(),
        mass,
        molar_mass,
        half_life,
        removal_time,
    )

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


def test_load_and_scale_mock() -> None:
    """Test that loading and scaling works on a mock spectrum."""
    # Create a fake spectrum
    energy = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3], dtype=float)
    dn_de = np.array([10, 20, 30, 40, 50, 60, 70], dtype=float)
    errors = np.array([1, 2, 3, 4, 5, 6, 7], dtype=float)
    data = np.column_stack((energy, dn_de, errors))

    # Define scaling parameters
    mass = 1.0  # kg
    molar_mass = 100.0  # g/mol
    half_life = 1.0  # year
    removal_time = 0.5  # years

    _test_load_and_scale(
        data,
        "test_isotope",
        mass,
        molar_mass,
        half_life,
        removal_time,
        max_energy=3,
    )


def _test_add_spec(spectra_list: ROOT.TList) -> ROOT.TH1D:
    """Test that adding multiple spectra works as expected."""
    # Add the spectra together
    combined_spec = add_spec(spectra_list)

    # Test basic properties
    assert isinstance(combined_spec, ROOT.TH1D), "Combined spectrum is not a TH1D"
    assert combined_spec.GetNbinsX() == spectra_list[0].GetNbinsX(), (
        "Combined spectrum has different number of bins"
    )

    # Test that each bin content and error has been added correctly
    for nbin in range(1, combined_spec.GetNbinsX() + 1):
        expected_content = sum(spec.GetBinContent(nbin) for spec in spectra_list)
        combined_content = combined_spec.GetBinContent(nbin)
        assert np.isclose(combined_content, expected_content), (
            f"Combined spectrum bin {nbin} content mismatch"
        )

        expected_error = np.sqrt(
            sum(spec.GetBinError(nbin) ** 2 for spec in spectra_list),
        )
        combined_error = combined_spec.GetBinError(nbin)
        assert np.isclose(combined_error, expected_error), (
            f"Combined spectrum bin {nbin} error mismatch"
        )

    return combined_spec


def test_add_spec_mock() -> None:
    """Test that adding spectra works on mock spectra."""
    # Create some fake spectra
    energy = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3], dtype=float)
    dn_de = np.array([10, 20, 30, 40, 50, 60, 70], dtype=float)
    errors = np.array([1, 2, 3, 4, 5, 6, 7], dtype=float)
    spec1 = create_spec(
        energy=energy,
        dn_de=dn_de,
        errors=errors,
        name="test_isotope1",
    )
    spec2 = create_spec(  # We'll make this one have double the counts and errors
        energy=energy,
        dn_de=dn_de * 2,
        errors=errors * 2,
        name="test_isotope2",
    )
    spectra_list = ROOT.TList()
    spectra_list.Add(spec1)
    spectra_list.Add(spec2)

    # Add the spectra together
    combined_spec = _test_add_spec(spectra_list)
    expected_contents = (
        dn_de + dn_de * 2
    )  # Sum of contents from both spectra, easy with numpy
    expected_errors = np.sqrt(errors**2 + (errors * 2) ** 2)  # Quadrature sum of errors
    for i, nbin in enumerate(range(1, combined_spec.GetNbinsX() + 1)):
        assert np.isclose(combined_spec.GetBinContent(nbin), expected_contents[i]), (
            f"Combined mock spectrum bin {nbin} content mismatch"
        )
        assert np.isclose(combined_spec.GetBinError(nbin), expected_errors[i]), (
            f"Combined mock spectrum bin {nbin} error mismatch"
        )


def test_sample_spec() -> None:
    """Test sample size and bounds from a ROOT spectrum."""
    ROOT.gRandom.SetSeed(1234)  # Fix seed for reproducibility

    # Create a fake spectrum
    spec = ROOT.TH1D("spec", "", 3, 0, 3)
    spec.SetBinContent(1, 1.0)
    spec.SetBinContent(2, 2.0)
    spec.SetBinContent(3, 3.0)

    # Test sampling 5 values from the spectrum
    samples = sample_spec(spec, samples=5)

    # Samples should be fixed given the seed
    samples_ref = [
        1.0745583504904062,
        1.9929909987840801,
        2.244217532686889,
        2.6356768854893744,
        1.81318321172148,
    ]
    assert samples == samples_ref, "Sampled values do not match reference"


def test_sample_spec_output() -> None:
    """Test CSV output writing when a filename is provided."""
    ROOT.gRandom.SetSeed(1234)  # Fix seed for reproducibility

    # Create a fake spectrum
    spec = ROOT.TH1D("spec", "", 3, 0, 3)
    spec.SetBinContent(1, 1.0)
    spec.SetBinContent(2, 2.0)
    spec.SetBinContent(3, 3.0)

    # Test sampling with output to a CSV file
    filename = Path("samples")
    samples = sample_spec(spec, samples=3, output_filename=filename)

    # Check that the file was created and contents are correct
    # (note that .csv extension is added automatically by sample_spec)
    csv_file = filename.with_suffix(".csv")
    saved = np.loadtxt(csv_file, delimiter=",")
    assert saved.tolist() == samples, "Saved samples do not match returned samples"

    # Clean up the file after testing
    csv_file.unlink()


def test_write_spec() -> None:
    """Test that spectra can be written to a ROOT file and read back correctly."""
    # Create a fake spectrum
    energy = np.array([0, 1, 2, 3, 4, 5, 6], dtype=float)
    dn_de = np.array([10, 20, 30, 40, 50, 60, 70], dtype=float)
    errors = np.array([1, 2, 3, 4, 5, 6, 7], dtype=float)
    spec = create_spec(
        energy=energy,
        dn_de=dn_de,
        errors=errors,
        name="test_isotope",
    )

    # Write the spectrum to a CSV
    filename = Path("test_spectrum")
    write_spec(spec, filename)
    energy, flux = _spec_to_arrays(spec)

    # Check that the file was created and contents are correct
    # (note that .csv extension is added automatically by sample_spec)
    csv_file = filename.with_suffix(".csv")
    data_loaded = np.loadtxt(
        csv_file,
        delimiter=",",
        skiprows=1,
    )
    energy_loaded = data_loaded[:, 0]
    flux_loaded = data_loaded[:, 1]
    assert np.array_equal(energy, energy_loaded), (
        "Loaded energy does not match original"
    )
    assert np.array_equal(flux, flux_loaded), "Loaded flux does not match original"

    # Clean up the file after testing
    csv_file.unlink()
