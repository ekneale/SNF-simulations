"""Define proportions of isotopes in SNF and calculate total antineutrino spectrum."""

from typing import NamedTuple

import numpy as np
import ROOT

from .load_data import load_antineutrino_data, load_isotope_data
from .spec import add_spec, load_spec

# TODO: these proportions could be loaded from a data file instead
SIZEWELL_PROPORTIONS = {
    "Sr90": 5.356e-4,
    "Y90": 1.3922e-7,
    "Pu241": 1.316e-3,
    "Cs137": 1.212e-3,
    "Am242": 3.554e-8,
    "Cs135": 3.1282e-4,
    "I129": 1.7535e-4,
    "Np239": 7.5852e-5,
    "Tc99": 7.9742e-4,
    "Zr93": 1.7681e-6,
    "Ce144": 4.0111e-4,
    "Kr88": 1.427e-10,
    "Pr144": 1.6896e-8,
    "Rb88": 1.6645e-11,
    "Rh106": 1.6389e-10,
    "Ru106": 1.7496e-4,
}
HARTLEPOOL_PROPORTIONS = {
    "Sr90": 4.2912e-4,
    "Y90": 1.0953e-7,
    "Pu241": 5.3075e-4,
    "Cs137": 8.6117e-4,
    "Am242": 1.3708e-8,
    "Cs135": 3.8379e-4,
    "I129": 1.2097e-4,
    "Np239": 2.0904e-5,
    "Tc99": 6.12e-4,
    "Zr93": 5.5068e-4,
    "Ce144": 1.7271e-4,
    "Kr88": 6.4374e-11,
    "Pr144": 7.2749e-9,
    "Rb88": 7.5089e-11,
    "Rh106": 5.7428e-11,
    "Ru106": 6.1306e-5,
}


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


def get_total_spec(
    cask_name: str,
    isotope_proportions: dict,
    total_mass: float = 1000,
    removal_time: float = 0,
    max_energy: float = 6000,
) -> ROOT.TH1D:
    """Calculate the total antineutrino spectrum from spent nuclear fuel.

    Args:
        cask_name: Name of the SNF cask.
        isotope_proportions: Dictionary of isotope proportions of the total mass.
        total_mass: Total mass of SNF (kg). Default is 1000 kg.
        removal_time: Time since removal from reactor (years). Default is 0.
        max_energy: Maximum energy to consider (keV). Default is 6000 keV.

    Returns:
        ROOT.TH1D: Total combined antineutrino spectrum as a ROOT histogram.

    """
    # Load the isotope data dicts
    isotopes = list(isotope_proportions.keys())
    molar_masses, half_lives = load_isotope_data(isotopes)
    isotope_data = load_antineutrino_data(isotopes)

    # Calculate the mass of each isotope from the input proportions of the total mass
    masses = {
        isotope: prop * total_mass for isotope, prop in isotope_proportions.items()
    }

    # Create the scaled spectra of each isotope and add them to a ROOT TList
    spectra = ROOT.TList()
    for isotope in isotopes:
        name = f"{isotope}{removal_time}{cask_name}"
        spec = load_spec(
            isotope_data[isotope],
            name,
            masses[isotope],
            molar_masses[isotope],
            half_lives[isotope],
            removal_time,
        )
        spectra.Add(spec)

    # Add any extra newly-created isotopes from decays.
    if removal_time != 0:
        # All of these decay chains have a branching ratio of 1.
        # If any additional isotopes were to be added with decay chains
        # involving more beta emitting isotopes then they can be added here.
        # TODO: work out how these are selected, if we can define them dynamically
        # or from an input file then that would be ideal.
        decay_chains = (
            DecayChain("Sr90", "Y90"),
            DecayChain("Ce144", "Pr144"),
            DecayChain("Kr88", "Rb88"),
            DecayChain("Ru106", "Rh106"),
        )

        for chain in decay_chains:
            if chain.parent not in masses or chain.daughter not in isotope_data:
                continue  # Skip if the isotope data is absent

            name = f"additional {chain.daughter}{removal_time}{cask_name}"
            daughter_mass = get_decay_mass(
                time_elapsed=removal_time,
                parent_mass=masses[chain.parent],
                parent_half_life=half_lives[chain.parent],
                daughter_half_life=half_lives[chain.daughter],
                branching_ratio=chain.branching_ratio,
            )
            spec = load_spec(
                isotope_data[chain.daughter],
                name,
                daughter_mass,
                molar_masses[chain.daughter],
                half_lives[chain.daughter],
                0,
            )
            spectra.Add(spec)

    # Sum all the spectra to get the total spectrum
    total_spec = add_spec(spectra)
    total_spec.SetTitle("Total Spectrum")
    total_spec.GetXaxis().SetRangeUser(0, max_energy)
    return total_spec
