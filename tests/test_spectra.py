"""Unit tests for loading antineutrino spectra data."""

import numpy as np
import ROOT
from snf_simulations.load_data import load_antineutrino_data
from snf_simulations.load_spec import load_spec


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
