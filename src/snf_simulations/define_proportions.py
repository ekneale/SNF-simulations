"""Define proportions of isotopes in SNF and calculate total antineutrino spectrum."""

from typing import NamedTuple

import numpy as np
import ROOT

from .load_data import load_antineutrino_data, load_isotope_data
from .spec import add_spec, load_spec


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


def TotSpec(
    cask_name,
    removal_time=0,
    total_mass=1000,
    max_E=6000,
    Sr90_prop=0,
    Y90_prop=0,
    Pu241_prop=0,
    Cs137_prop=0,
    Am242_prop=0,
    Cs135_prop=0,
    I129_prop=0,
    Np239_prop=0,
    Tc99_prop=0,
    Zr93_prop=0,
    Ce144_prop=0,
    Kr88_prop=0,
    Pr144_prop=0,
    Rb88_prop=0,
    Rh106_prop=0,
    Ru106_prop=0,
):
    # total mass is in kg
    # all times given in years

    # inputted proportions of isotopes defined by the user
    proportions = {
        "Sr90": Sr90_prop,
        "Y90": Y90_prop,
        "Pu241": Pu241_prop,
        "Cs137": Cs137_prop,
        "Am242": Am242_prop,
        "Cs135": Cs135_prop,
        "I129": I129_prop,
        "Np239": Np239_prop,
        "Tc99": Tc99_prop,
        "Zr93": Zr93_prop,
        "Ce144": Ce144_prop,
        "Kr88": Kr88_prop,
        "Pr144": Pr144_prop,
        "Rb88": Rb88_prop,
        "Rh106": Rh106_prop,
        "Ru106": Ru106_prop,
    }

    # Load the isotope data dicts
    isotopes = list(proportions.keys())
    molar_masses, half_lives = load_isotope_data(isotopes)
    isotope_data = load_antineutrino_data(isotopes)

    # Calculate the mass of each isotope from the input proportions of the total mass
    masses = {isotope: prop * total_mass for isotope, prop in proportions.items()}

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
    total_spec.GetXaxis().SetRangeUser(0, max_E)
    return total_spec
