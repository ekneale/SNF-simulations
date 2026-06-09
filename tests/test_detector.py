"""Unit tests for the Detector class."""

import numpy as np
import pytest

from snf_simulations.detector import Detector
from snf_simulations.physics import calculate_flux_at_distance
from snf_simulations.spec import Spectrum

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers


def test_create_detector() -> None:
    """Test basic Detector construction with valid values."""
    detector = Detector(volume=1.2, proton_density=4.6e22, name="test_detector")

    assert detector.volume == 1.2
    assert detector.proton_density == 4.6e22
    assert detector.name == "test_detector"
    assert np.isclose(detector.number_of_protons, 1.2 * 1e6 * 4.6e22)


def test_init() -> None:
    """Test Detector input validation."""
    with pytest.raises(ValueError, match="Detector volume must be a positive value"):
        Detector(volume=0.0, proton_density=4.6e22)

    with pytest.raises(ValueError, match="Detector volume must be a positive value"):
        Detector(volume=-1.0, proton_density=4.6e22)

    with pytest.raises(ValueError, match="Proton density must be a positive value"):
        Detector(volume=1.0, proton_density=0.0)

    with pytest.raises(ValueError, match="Proton density must be a positive value"):
        Detector(volume=1.0, proton_density=-1.0)


def test_repr() -> None:
    """Test repr for initialized and uninitialized Detector objects."""
    detector = Detector(volume=1.2, proton_density=4.6e22, name="test_detector")

    expected_repr = (
        f'<Detector "{detector.name}": '
        f"volume={detector.volume:.3e} m3, "
        f"proton_density={detector.proton_density:.3e} protons/cm3>"
    )
    assert repr(detector) == expected_repr

    uninitialized = Detector.__new__(Detector)
    assert repr(uninitialized) == "<Detector (uninitialized)>"

    unnamed = Detector(volume=1.2, proton_density=4.6e22)
    expected_repr_unnamed = (
        f"<Detector: "
        f"volume={unnamed.volume:.3e} m3, "
        f"proton_density={unnamed.proton_density:.3e} protons/cm3>"
    )
    assert repr(unnamed) == expected_repr_unnamed


def test_calculate_event_rate() -> None:
    """Test event rate calculation."""
    detector_volume = 2 * 0.6  # m3
    proton_density = 4.6e22  # protons/cm3
    detector = Detector(volume=detector_volume, proton_density=proton_density)

    # Use mock spectrum
    energy = np.array([0.0, 1000.0, 2000.0, 3000.0])
    flux = np.array([10.0, 20.0, 30.0])
    errors = np.zeros_like(flux)
    spec = Spectrum(energy=energy, flux=flux, errors=errors)
    event_rate = detector.calculate_event_rate(spec=spec, distance=40.0)

    number_of_protons = detector_volume * 1e6 * proton_density
    flux = calculate_flux_at_distance(spec.integrate(1806), distance=40.0)
    expected_event_rate = number_of_protons * 1e-44 * flux

    assert np.isclose(event_rate, expected_event_rate), (
        f"Calculated event rate ({event_rate}) does not match expected value "
        f"({expected_event_rate})"
    )

    # Test the efficiency parameter
    event_rate_half = detector.calculate_event_rate(
        spec=spec, distance=40.0, efficiency=0.5
    )
    assert np.isclose(event_rate_half, expected_event_rate * 0.5), (
        "Event rate with 50% efficiency should be half of the full efficiency rate"
    )

    # Test with a spectrum that has no flux above the IBD threshold (1806 keV)
    energy = np.array([0.0, 100.0, 200.0, 300.0])
    flux = np.array([10.0, 20.0, 30.0])
    errors = np.zeros_like(flux)
    spec = Spectrum(energy=energy, flux=flux, errors=errors)
    with pytest.raises(
        ValueError, match="Spectrum energy range does not cover the IBD threshold"
    ):
        detector.calculate_event_rate(spec=spec, distance=40.0)
