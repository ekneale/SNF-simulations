"""Calculate antineutrino spectra for spent nuclear fuel casks."""

from copy import deepcopy

import ROOT

from .data import load_antineutrino_data, load_isotope_data
from .physics import DecayChain, get_decay_mass, get_isotope_activity
from .spec import Spectrum, add_spec, load_spec


class Cask:
    """Class representing a cask of spent nuclear fuel.

    Attributes:
        isotope_proportions: The proportions of each isotope in the cask.
            Should be a dictionary where keys are isotope names and values are the
            proportion of the total mass that isotope represents (between 0 and 1).
        total_mass: The total mass of the cask in kg.
        name: The name of the cask.

    """

    def __init__(
        self,
        isotope_proportions: dict[str, float],
        total_mass: float = 1000,
        name: str = "Cask",
    ) -> None:
        """Initialize the Cask object."""
        self.name = name
        self.isotope_proportions = isotope_proportions
        self.total_mass = total_mass

        if not self.isotope_proportions:
            msg = "isotope_proportions must not be empty"
            raise ValueError(msg)
        if any(proportion < 0 for proportion in self.isotope_proportions.values()):
            msg = "isotope_proportions values must be non-negative"
            raise ValueError(msg)
        if self.total_mass <= 0:
            msg = "total_mass must be positive"
            raise ValueError(msg)

        # Store all constant isotope data
        self.isotopes = list(self.isotope_proportions.keys())
        self.isotope_masses = {
            isotope: proportion * self.total_mass
            for isotope, proportion in self.isotope_proportions.items()
        }
        # TODO: use medvedev package for molar masses and half-lives
        self.molar_masses, self.half_lives = load_isotope_data(self.isotopes)
        self.isotope_spectra = {
            isotope: Spectrum.from_isotope(isotope) for isotope in self.isotopes
        }

    def __repr__(self) -> str:
        """Return a string representation of the Cask object."""
        try:
            repr_str = f'<Cask "{self.name}", total_mass={self.total_mass} kg>'
        except AttributeError:
            return "<Cask (uninitialized)>"
        else:
            return repr_str

    def _get_component_spectra(self, removal_time: float = 0) -> list[Spectrum]:
        """Get the individual antineutrino spectra for each isotope in the cask.

        Returns:
            A list of Spectrum objects, representing the antineutrino spectra for each
            isotope as well as any additional isotopes created from decays since
            removal from the reactor.

        """
        if removal_time < 0:
            msg = "removal_time must be non-negative"
            raise ValueError(msg)

        # Get the antineutrino spectra for each isotope
        spectra = []
        for isotope in self.isotopes:
            # Get the antineutrino spectrum
            spec = self.isotope_spectra[isotope]

            # Scale based on given removal time
            activity = get_isotope_activity(
                mass=self.isotope_masses[isotope],
                molar_mass=self.molar_masses[isotope],
                half_life=self.half_lives[isotope],
                removal_time=removal_time,
            )
            scaled_spec = spec * activity
            spectra.append(scaled_spec)

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
                if (
                    chain.parent not in self.isotopes
                    or self.isotope_masses[chain.parent] == 0
                ):
                    # No parent isotope in the cask, so skip this decay
                    continue

                if chain.daughter not in self.isotopes:
                    # We won't have the spectrum or hl/mm data cached
                    daughter_spec = Spectrum.from_isotope(chain.daughter)
                    # TODO: load_isotope_data should take a single isotope name
                    _molar_masses, _half_lives = load_isotope_data([chain.daughter])
                    daughter_molar_mass = _molar_masses[chain.daughter]
                    daughter_half_life = _half_lives[chain.daughter]
                else:
                    daughter_spec = deepcopy(self.isotope_spectra[chain.daughter])
                    daughter_molar_mass = self.molar_masses[chain.daughter]
                    daughter_half_life = self.half_lives[chain.daughter]
                daughter_spec.name = f"{chain.parent}->{chain.daughter}"

                # Calculate the mass of the daughter isotope
                daughter_mass = get_decay_mass(
                    time_elapsed=removal_time,
                    parent_mass=self.isotope_masses[chain.parent],
                    parent_half_life=self.half_lives[chain.parent],
                    daughter_half_life=daughter_half_life,
                    branching_ratio=chain.branching_ratio,
                )

                # Scale based on the daughter mass
                # Note the removal time is set to 0 here, since the daughter mass
                # already accounts for decay during the time since the cask was removed
                activity = get_isotope_activity(
                    mass=daughter_mass,
                    molar_mass=daughter_molar_mass,
                    half_life=daughter_half_life,
                    removal_time=0,
                )
                scaled_spec = daughter_spec * activity
                spectra.append(scaled_spec)
        return spectra

    def get_total_spectrum(self, removal_time: float = 0) -> Spectrum:
        """Calculate the total antineutrino spectrum as a Spectrum object.

        Args:
            removal_time: The time in years since the cask was removed from the reactor.

        """
        spectra = self._get_component_spectra(removal_time)

        # Equalise all the spectra to 1keV bins to allow combining,
        # going from 0 to the maximum energy across all spectra.
        max_energy = max(spec.energy[-1] for spec in spectra)
        for spec in spectra:
            spec.equalise(width=1, min_energy=0, max_energy=max_energy)

        # Sum all the spectra to get the total cask spectrum
        total_spec = spectra[0]
        for spec in spectra[1:]:
            total_spec = total_spec + spec
        total_spec.name = f"{self.name} total spectrum"
        return total_spec


def get_total_spec(
    cask_name: str,
    isotope_proportions: dict,
    total_mass: float = 1000,
    removal_time: float = 0,
    max_energy: float | None = None,
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
            max_energy=max_energy,
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
            if chain.parent not in masses:
                continue  # Skip if the isotope data is absent

            if chain.daughter not in isotopes:
                daughter_data = load_antineutrino_data([chain.daughter])
                daughter_data = daughter_data[chain.daughter]
                # We won't have the spectrum or hl/mm data cached
                _molar_masses, _half_lives = load_isotope_data([chain.daughter])
                daughter_molar_mass = _molar_masses[chain.daughter]
                daughter_half_life = _half_lives[chain.daughter]
            else:
                daughter_data = isotope_data[chain.daughter]
                daughter_molar_mass = molar_masses[chain.daughter]
                daughter_half_life = half_lives[chain.daughter]

            name = f"additional {chain.daughter}{removal_time}{cask_name}"
            daughter_mass = get_decay_mass(
                time_elapsed=removal_time,
                parent_mass=masses[chain.parent],
                parent_half_life=half_lives[chain.parent],
                daughter_half_life=daughter_half_life,
                branching_ratio=chain.branching_ratio,
            )
            spec = load_spec(
                daughter_data,
                name,
                daughter_mass,
                daughter_molar_mass,
                daughter_half_life,
                0,
                max_energy=max_energy,
            )
            spectra.Add(spec)

    # Sum all the spectra to get the total spectrum
    total_spec = add_spec(spectra)
    total_spec.SetTitle("Total Spectrum")
    if max_energy is not None:
        total_spec.GetXaxis().SetRangeUser(0, max_energy)
    return total_spec
