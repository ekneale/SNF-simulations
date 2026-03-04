"""Unit tests for the Cask class."""

import pytest

from snf_simulations.cask import Cask
from snf_simulations.spec import Spectrum

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts


def test_create_cask() -> None:
    """Test basic Cask construction with default values."""
    isotope_proportions = {"Sr90": 0.5, "Cs137": 0.5}
    total_mass = 2000
    name = "test_cask"

    # Create the Cask
    cask = Cask(isotope_proportions, total_mass, name=name)

    # Test basic properties
    assert cask.isotope_proportions == isotope_proportions
    assert cask.total_mass == total_mass
    assert cask.name == name

    # Test derived properties
    assert cask.isotopes == ["Sr90", "Cs137"]
    assert cask.isotope_masses == {"Sr90": 1000.0, "Cs137": 1000.0}


def test_init() -> None:
    """Test that Cask constructor validates inputs."""
    with pytest.raises(ValueError, match="isotope_proportions must not be empty"):
        Cask(isotope_proportions={}, total_mass=1000)

    with pytest.raises(
        ValueError, match="isotope_proportions values must be non-negative"
    ):
        Cask(isotope_proportions={"Sr90": -1.0}, total_mass=1000)

    with pytest.raises(ValueError, match="total_mass must be positive"):
        Cask(isotope_proportions={"Sr90": 1.0}, total_mass=-1000)

    with pytest.raises(ValueError, match="total_mass must be positive"):
        Cask(isotope_proportions={"Sr90": 1.0}, total_mass=0)


def test_repr() -> None:
    """Test repr for initialized and uninitialized Cask objects."""
    isotope_proportions = {"Sr90": 0.5, "Cs137": 0.5}
    total_mass = 2000
    name = "test_cask"
    cask = Cask(isotope_proportions, total_mass, name=name)

    expected_repr = f'<Cask "{cask.name}", total_mass={cask.total_mass} kg>'
    assert repr(cask) == expected_repr, (
        "Cask repr should include cask name and total mass"
    )

    # Create an instance without calling __init__ to exercise fallback repr path.
    uninitialized = Cask.__new__(Cask)
    assert repr(uninitialized) == "<Cask (uninitialized)>", (
        "Cask repr should handle uninitialized objects gracefully"
    )


def test_get_total_spectrum() -> None:
    """Test that get_total_spectrum returns a Spectrum object."""
    isotope_proportions = {"Sr90": 0.5, "Cs137": 0.5}
    total_mass = 2000
    name = "test_cask"
    cask = Cask(isotope_proportions, total_mass, name=name)

    # We won't test the actual spectrum content here, just that it returns a Spectrum
    # object with the expected name.
    spec = cask.get_total_spectrum(removal_time=1.0)
    assert isinstance(spec, Spectrum), "get_total_spectrum should return a Spectrum"
    assert spec.name == f"{cask.name} total spectrum", (
        "Spectrum name should be based on cask name"
    )


def test_get_total_spectrum_inputs() -> None:
    """Test that get_total_spectrum validates its inputs."""
    isotope_proportions = {"Sr90": 0.5, "Cs137": 0.5}
    total_mass = 2000
    name = "test_cask"
    cask = Cask(isotope_proportions, total_mass, name=name)

    with pytest.raises(ValueError, match="removal_time must be non-negative"):
        cask.get_total_spectrum(removal_time=-1.0)
