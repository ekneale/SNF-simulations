"""Unit tests for the Cask class."""

from pathlib import Path

import pytest

from snf_simulations.cask import Cask
from snf_simulations.data import get_isotope_properties
from snf_simulations.spec import Spectrum

from .test_data_fispin import _write_tabqfile

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts


def test_create_cask() -> None:
    """Test basic Cask construction with default values."""
    isotope_masses = {"Sr90": 1000.0, "Cs137": 1000.0}  # kg
    initial_cooling_time = 10.0  # years
    name = "test_cask"

    # Create the Cask
    cask = Cask(isotope_masses, initial_cooling_time, name)

    # Test basic properties
    assert cask.isotope_masses == isotope_masses
    assert cask.name == name
    assert cask.initial_cooling_time == initial_cooling_time

    # Test derived properties
    assert cask.isotopes == ["Sr90", "Cs137"]
    assert cask.isotope_properties == {
        isotope: get_isotope_properties(isotope) for isotope in cask.isotopes
    }
    assert cask.isotope_spectra == {
        isotope: Spectrum.from_isotope(isotope) for isotope in cask.isotopes
    }


def test_init() -> None:
    """Test that Cask constructor validates inputs."""
    with pytest.raises(ValueError, match="isotope_masses must not be empty"):
        Cask(isotope_masses={}, name="test_cask")

    with pytest.raises(ValueError, match="isotope_masses values must be non-negative"):
        Cask(isotope_masses={"Sr90": -1.0}, name="test_cask")

    with pytest.raises(ValueError, match="initial_cooling_time must be non-negative"):
        Cask(isotope_masses={"Sr90": 1000.0}, initial_cooling_time=-1.0)


def test_repr() -> None:
    """Test repr for initialized and uninitialized Cask objects."""
    isotope_masses = {"Sr90": 1000.0, "Cs137": 1000.0}  # kg
    initial_cooling_time = 10.0  # years
    name = "test_cask"
    cask = Cask(isotope_masses, initial_cooling_time, name)

    expected_repr = (
        f'<Cask "{name}": '
        f"{len(isotope_masses)} isotopes, "
        f"cooling time={initial_cooling_time:.3e} years>"
    )
    assert repr(cask) == expected_repr, (
        "Cask repr should include cask name, total mass and cooling time"
    )

    # Create an instance without calling __init__ to exercise fallback repr path.
    uninitialized = Cask.__new__(Cask)
    assert repr(uninitialized) == "<Cask (uninitialized)>", (
        "Cask repr should handle uninitialized objects gracefully"
    )

    # Create an unnamed instance.
    unnamed = Cask(isotope_masses, initial_cooling_time)
    expected_repr = (
        f"<Cask: "
        f"{len(isotope_masses)} isotopes, "
        f"cooling time={initial_cooling_time:.3e} years>"
    )
    assert repr(unnamed) == expected_repr, (
        "Cask repr should include total mass and cooling time "
        "even if name is not provided"
    )


def test_from_tabqfile(tmp_path: Path) -> None:
    """Test that Cask.from_tabqfile creates a Cask with expected properties."""
    filepath = _write_tabqfile(tmp_path)

    cask = Cask.from_tabqfile(filepath, total_mass=1000.0, name="test_cask")

    assert isinstance(cask, Cask), "from_tabqfile should return a Cask instance"
    assert cask.name == "test_cask", "Cask name should match the provided name"
    assert cask.isotopes == ["Sr90", "Cs137"], (
        "Isotopes should be extracted from the .tbQ file and capitalized"
    )
    # It should have selected the lowest cooling time by default (12 hours)
    assert cask.initial_cooling_time == pytest.approx(0.5 / 365.2425)
    expected_isotope_masses = {"Sr90": 500.0, "Cs137": 500.0}
    assert cask.isotope_masses == expected_isotope_masses, "Incorrect isotope masses"


def test_from_tabqfile_default_name(tmp_path: Path) -> None:
    """Test that from_tabqfile can select a subset of isotopes."""
    filepath = _write_tabqfile(tmp_path)

    cask = Cask.from_tabqfile(filepath, total_mass=1000.0)
    assert cask.name == "sample", "Cask name should be set by the file name"

    cask = Cask.from_tabqfile(str(filepath), total_mass=1000.0)
    assert cask.name == "sample", "Cask name should be set by the file name"


def test_from_tabqfile_selected_isotopes(tmp_path: Path) -> None:
    """Test that from_tabqfile can select a subset of isotopes."""
    filepath = _write_tabqfile(tmp_path)

    cask = Cask.from_tabqfile(
        filepath, total_mass=1000.0, name="test_cask", isotopes=["Sr90"]
    )

    assert cask.isotopes == ["Sr90"], "Cask should only include the selected isotope"
    assert cask.isotope_masses == {"Sr90": 500.0}, "Incorrect isotope masses"


def test_get_total_spectrum() -> None:
    """Test that get_total_spectrum returns a Spectrum object."""
    isotope_masses = {"Sr90": 1000.0, "Cs137": 1000.0}  # kg
    initial_cooling_time = 10.0  # years
    name = "test_cask"
    cask = Cask(isotope_masses, initial_cooling_time, name=name)

    # We won't test the actual spectrum content here, just that it returns a Spectrum
    # object with the expected name.
    spec = cask.get_total_spectrum(cooling_time=20.0)
    assert isinstance(spec, Spectrum), "get_total_spectrum should return a Spectrum"
    assert spec.name == f"{cask.name}", (
        "Spectrum name should be the same as the Cask name"
    )


def test_get_total_spectrum_inputs() -> None:
    """Test that get_total_spectrum validates its inputs."""
    isotope_masses = {"Sr90": 1000.0, "Cs137": 1000.0}  # kg
    initial_cooling_time = 10.0  # years
    name = "test_cask"
    cask = Cask(isotope_masses, initial_cooling_time, name=name)

    with pytest.raises(ValueError, match="cooling_time must be non-negative"):
        cask.get_total_spectrum(cooling_time=-1.0)
    with pytest.raises(
        ValueError, match="cannot be less than the initial cask cooling time"
    ):
        cask.get_total_spectrum(cooling_time=initial_cooling_time - 1.0)


def test_get_total_spectrum_daughter() -> None:
    """Test that get_total_spectrum includes daughter isotope spectra."""
    isotope_masses = {"Sr90": 1000.0, "Y90": 1000.0}  # kg
    initial_cooling_time = 10.0  # years
    name = "test_cask"
    cask = Cask(isotope_masses, initial_cooling_time, name=name)

    component_spec = cask._get_component_spectra(cooling_time=20.0)
    assert len(component_spec) == len(isotope_masses) + 1, (
        "Should have spectra for both input isotopes and the SR90->Y90 daughter"
    )
    component_spec_names = [spec.name for spec in component_spec]
    assert "Sr90" in component_spec_names, "Should have a spectrum for Sr90"
    assert "Y90" in component_spec_names, "Should have a spectrum for Y90"
    assert "Sr90->Y90" in component_spec_names, (
        "Should have a spectrum for the Sr90 -> Y90 decay"
    )
