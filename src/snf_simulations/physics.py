"""Module for atomic physics calculations."""

from typing import NamedTuple

import numpy as np


class DecayChain(NamedTuple):
    """Class to represent a decay chain from parent to daughter isotope."""

    parent: str
    daughter: str
    branching_ratio: float = 1.0


def get_isotope_activity(
    time_elapsed: float,
    mass: float,
    molar_mass: float,
    half_life: float,
) -> float:
    """Calculate the activity of an isotope after a given time.

    Args:
        time_elapsed: Time elapsed in years.
        mass: Mass of the isotope in kg.
        molar_mass: Molar mass of the isotope in g/mol.
        half_life: Half-life of the isotope in years.

    Returns:
        Activity of the isotope in decays per second (Becquerels)
        after the given time.

    """
    # Convert mass to number of atoms (kg to g, then to moles, then to atoms)
    number_of_atoms = (mass * 1000 / molar_mass) * 6.022e23
    # Calculate initial activity (in decays per second aka Becquerels)
    lambda_ = np.log(2) / (half_life * 365 * 24 * 60 * 60)
    initial_activity = number_of_atoms * lambda_
    # Calculate activity after time (still in decays per second)
    activity = initial_activity * np.exp(
        -1 * lambda_ * time_elapsed * 365 * 24 * 60 * 60,
    )
    return activity


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
        time_elapsed: Time elapsed since initial measurement in years.
        parent_mass: Initial mass of parent isotope in kg.
        parent_half_life: Half-life of parent isotope in years.
        daughter_half_life: Half-life of daughter isotope in years.
        branching_ratio: Branching ratio for this decay chain.

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


def calculate_flux_at_distance(
    total_flux: float,
    distance: float,
) -> float:
    """Calculate the antineutrino flux at a given distance from the source.

    We assume the source is a point source emitting isotropically in all directions,
    so the flux decreases with the square of the distance from the source.

    Args:
        total_flux: The total antineutrino flux emitted by the source in s^-1 .
        distance: Distance from the source in meters.

    Returns:
        Antineutrino flux at the given distance in cm^-2 s^-1.

    """
    # Note we convert distance to cm
    flux = (1 / (4 * np.pi * (distance * 100) ** 2)) * (total_flux)
    return flux
