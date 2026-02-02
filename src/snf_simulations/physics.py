"""Module for atomic physics calculations."""

from typing import NamedTuple

import numpy as np


class DecayChain(NamedTuple):
    """Class to represent a decay chain from parent to daughter isotope."""

    parent: str
    daughter: str
    branching_ratio: float = 1.0


def get_decay_mass(
    time_elapsed: float,
    parent_mass: float,
    parent_half_life: float,
    daughter_half_life: float,
    branching_ratio: float = 1,
) -> float:
    """Calculate the mass of a daughter isotope created from parent decay.

    Computes the mass of a daughter isotope that has been created from the
    radioactive decay of its parent isotope using first-order decay equations.

    Args:
        time_elapsed: Time elapsed since initial measurement (years).
        parent_mass: Initial mass of parent isotope (kg).
        parent_half_life: Half-life of parent isotope (years).
        daughter_half_life: Half-life of daughter isotope (years).
        branching_ratio: Branching ratio for this decay chainway. Defaults to 1.

    Returns:
        Mass of the daughter isotope (kg).

    """
    # Decay constants (natural log of 2 divided by half-life)
    parent_decay_constant = np.log(2) / parent_half_life
    daughter_decay_constant = np.log(2) / daughter_half_life

    # Bateman equation for daughter isotope mass
    daughter_mass = (
        branching_ratio
        * (parent_decay_constant / (daughter_decay_constant - parent_decay_constant))
        * parent_mass
        * (
            np.exp(-parent_decay_constant * time_elapsed)
            - np.exp(-daughter_decay_constant * time_elapsed)
        )
    )
    return daughter_mass


def calculate_flux(spec, distance):
    """Calculate the antineutrino flux at a given distance from the source.

    Args:
        spec (ROOT.TH1D): The total antineutrino spectrum histogram.
        distance (float): Distance from the source in meters.

    Returns:
        flux (float): Antineutrino flux at the given distance in cm^{-2} s^{-1}.

    """
    # The IBD reaction threshold is 1800 keV, so we integrate above this energy.
    # This gives the total flux of antineutrinos per second above the threshold.
    total_flux = spec.Integral(1801, 6000)

    # Calculate the flux at the given distance, assuming it's a point source emitting
    # isotropically in all directions.
    # Note we convert to cm, to get the flux is in cm^-2 s^-1
    flux = (1 / (4 * np.pi * (distance * 100) ** 2)) * (total_flux)
    return flux


def calculate_event_rate(
    flux,
    lower_efficiency=0.2,
    upper_efficiency=0.4,
):
    """Calculate the expected event rate in the VIDARR detector for a given flux.

    Args:
        flux (float): Antineutrino flux in cm^-2 s^-1.
        lower_efficiency (float): Lower limit on detection efficiency. Default is 0.2.
        upper_efficiency (float): Upper limit on detection efficiency. Default is 0.4

    Returns:
        rate_lower, rate_upper (tuple of floats): event rates in s^{-1}
            for the lower and upper efficiency limits.

    """
    # Calculate the number of target protons in the detector
    # VIDARR is a plastic scintillator detector with a volume of (1.52 x 1.52 x 0.7) m^3
    detector_volume = 1.52 * 1.52 * 0.7  # m^3
    detector_volume = detector_volume * 1e6  # convert to cm^3
    proton_density = 4.6e22  # number density of protons in cm^-3
    number_of_protons = detector_volume * proton_density

    # Calculate the event rate using the flux and number of target protons
    cross_section = 1e-44  # cm^2, approximate IBD cross-section
    event_rate = number_of_protons * cross_section * flux

    # Apply detection efficiency
    rate_lower = event_rate * lower_efficiency
    rate_upper = event_rate * upper_efficiency
    return rate_lower, rate_upper
