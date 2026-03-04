"""Calculate antineutrino spectra for spent nuclear fuel casks."""

import ROOT

from .data import load_antineutrino_data, load_isotope_data
from .physics import DecayChain, get_decay_mass
from .spec import add_spec, load_spec


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
        total_mass: Total mass of SNF (kg).
        removal_time: Time since removal from reactor (years).
        max_energy: Maximum energy to consider (keV).

    Returns:
        Total combined antineutrino spectrum as a ROOT histogram.

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
