"""Calculate antineutrino spectra for spent nuclear fuel casks."""

from collections.abc import Collection
from copy import deepcopy
from pathlib import Path
from typing import cast

from .data import get_isotope_masses, get_isotope_properties
from .physics import DecayChain, get_decay_mass, get_isotope_activity
from .spec import Spectrum

# Define a default list of isotopes to include in the cask spectrum if the user doesn't
# specify their own list.
DEFAULT_ISOTOPES = [
    "Sr90",
    "Y90",
    "Pu241",
    "Cs137",
    "Am242",
    "Cs135",
    "I129",
    "Np239",
    "Tc99",
    "Zr93",
    "Ce144",
    "Kr88",
    "Pr144",
    "Rb88",
    "Rh106",
    "Ru106",
]
_DEFAULT_ISOTOPES = object()  # Sentinel value (see Cask.from_tabqfile)


class Cask:
    """Class representing a cask of spent nuclear fuel.

    Attributes:
        isotope_masses: The masses of each isotope in the cask.
            Should be a dictionary where keys are isotope names and values are the
            mass of the isotope in the cask (in kg).
        initial_cooling_time: The time since the cask was removed from the reactor,
            in years, corresponding to the time that the isotope masses were calculated.
            This initial age will be subtracted from the requested cooling times when
            calculating the antineutrino spectrum, to account for decay during the time
            since removal.
        name: An optional name for the cask.

    """

    def __init__(
        self,
        isotope_masses: dict[str, float],
        initial_cooling_time: float = 0,
        name: str | None = None,
    ) -> None:
        """Initialize the Cask object."""
        self.isotope_masses = isotope_masses
        self.initial_cooling_time = initial_cooling_time
        self.name = name

        if not self.isotope_masses:
            msg = "isotope_masses must not be empty"
            raise ValueError(msg)
        if any(mass < 0 for mass in self.isotope_masses.values()):
            msg = "isotope_masses values must be non-negative"
            raise ValueError(msg)
        if self.initial_cooling_time < 0:
            msg = "initial_cooling_time must be non-negative"
            raise ValueError(msg)

        # Store all constant isotope data
        self.isotopes = list(self.isotope_masses.keys())
        self.isotope_properties = {
            isotope: get_isotope_properties(isotope) for isotope in self.isotopes
        }
        self.isotope_spectra = {
            isotope: Spectrum.from_isotope(isotope) for isotope in self.isotopes
        }

    def __repr__(self) -> str:
        """Return a string representation of the Cask object."""
        try:
            name_str = f' "{self.name}"' if self.name is not None else ""
            repr_str = (
                f"<Cask{name_str}: "
                f"{len(self.isotope_masses)} isotopes, "
                f"cooling time={self.initial_cooling_time:.3e} years>"
            )
        except AttributeError:
            return "<Cask (uninitialized)>"
        else:
            return repr_str

    @classmethod
    def from_tabqfile(
        cls,
        filepath: str | Path,
        total_mass: float | None = None,
        isotopes: Collection[str] | None | object = _DEFAULT_ISOTOPES,
        time_str: str | None = None,
        name: str | None = None,
    ) -> "Cask":
        """Create a Cask object from a FISPIN .tbQ output file.

        Args:
            filepath: Path to the file to load.
            total_mass: The total mass of the cask to simulate (in kg).
                If None, the mass from the simulation file is used.
                If given, the isotopes will be scaled in proportion to the
                simulation mass.
            isotopes: Optional list of isotopes to include from the file.
                Defaults to a list of selected isotopes (see DEFAULT_ISOTOPES).
                If None, includes all isotopes in the file.
            time_str: Specific simulation time to extract data for.
                If None, the smallest time in the file is used.
                (see data.get_isotope_masses for details)
            name: Optional name for the cask.
                If None, a name is generated from the filename.

        Returns:
            A Cask object with the isotope masses loaded from the file.

        """
        isotope_masses, cooling_time = get_isotope_masses(filepath, time_str)

        # If a specific mass is given, scale the isotopes proportionally.
        if total_mass is not None:
            scaling_factor = total_mass / sum(isotope_masses.values())
            isotope_masses = {
                isotope: mass * scaling_factor
                for isotope, mass in isotope_masses.items()
            }

        # Filter isotopes
        # We use a sentinel value for the default here to allow distinguishing between
        # the user explicitly passing None (to include all isotopes) and
        # not passing anything (to use the default list of selected isotopes).
        # The cast is needed to satisfy type checking.
        # Note PEP 661 proposes a builtin sentinel class which would be helpful here.
        if isotopes is not None:
            if isotopes is _DEFAULT_ISOTOPES:
                selected_isotopes = DEFAULT_ISOTOPES
            else:
                selected_isotopes = cast(Collection[str], isotopes)
            isotope_masses = {
                isotope: mass
                for isotope, mass in isotope_masses.items()
                if isotope in selected_isotopes
            }
        if name is None:
            if isinstance(filepath, Path):
                name = filepath.stem
            elif filepath.endswith(".tbQ"):
                name = filepath.rsplit("/", maxsplit=1)[-1].split(".", maxsplit=1)[0]
        return cls(
            isotope_masses=isotope_masses,
            initial_cooling_time=cooling_time,
            name=name,
        )

    def _get_component_spectra(self, cooling_time: float = 0) -> list[Spectrum]:
        """Get the individual antineutrino spectra for each isotope in the cask.

        Args:
            cooling_time: The time in years since the cask was removed from the reactor.
                Note that this has to be greater than or equal to the
                initial_cooling_time of the cask.

        Returns:
            A list of Spectrum objects, representing the antineutrino spectra for each
            isotope as well as any additional isotopes created from decays since
            removal from the reactor.

        """
        if cooling_time < 0:
            msg = "cooling_time must be non-negative"
            raise ValueError(msg)
        if cooling_time < self.initial_cooling_time:
            msg = f"cooling_time ({cooling_time:.3e}) cannot be less than "
            msg += f"the initial cask cooling time ({self.initial_cooling_time:.3e})"
            raise ValueError(msg)

        # Take off the initial age of the cask, as the isotope masses should already
        # account for some initial decay since removal from the core.
        time_elapsed = cooling_time - self.initial_cooling_time

        # Get the antineutrino spectra for each isotope
        spectra = []
        for isotope in self.isotopes:
            # Get the antineutrino spectrum for this isotope
            spec = self.isotope_spectra[isotope]

            # Scale based on the time since the initial removal
            activity = get_isotope_activity(
                time_elapsed=time_elapsed,
                mass=self.isotope_masses[isotope],
                molar_mass=self.isotope_properties[isotope]["molar_mass"],
                half_life=self.isotope_properties[isotope]["half_life"],
            )
            spec *= activity
            spectra.append(spec)

        # Add any extra newly-created isotopes from decays since
        # the initial cooling time.
        if time_elapsed > 0:
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

                if chain.daughter not in self.isotope_properties:
                    # We won't have the spectrum or hl/mm data cached
                    daughter_spec = Spectrum.from_isotope(chain.daughter)
                    daughter_properties = get_isotope_properties(chain.daughter)
                else:
                    daughter_spec = deepcopy(self.isotope_spectra[chain.daughter])
                    daughter_properties = self.isotope_properties[chain.daughter]
                daughter_molar_mass = daughter_properties["molar_mass"]
                daughter_half_life = daughter_properties["half_life"]
                daughter_spec.name = f"{chain.parent}->{chain.daughter}"

                # Calculate the mass of the daughter isotope
                daughter_mass = get_decay_mass(
                    time_elapsed=time_elapsed,
                    parent_mass=self.isotope_masses[chain.parent],
                    parent_half_life=self.isotope_properties[chain.parent]["half_life"],
                    daughter_half_life=daughter_half_life,
                    branching_ratio=chain.branching_ratio,
                )

                # Scale the spectrum based on the daughter mass
                # Note the time_elapsed is set to 0 here, since the daughter mass
                # already accounts for any decay.
                activity = get_isotope_activity(
                    time_elapsed=0,
                    mass=daughter_mass,
                    molar_mass=daughter_molar_mass,
                    half_life=daughter_half_life,
                )
                scaled_spec = daughter_spec * activity
                spectra.append(scaled_spec)
        return spectra

    def get_total_spectrum(self, cooling_time: float = 0) -> Spectrum:
        """Calculate the total antineutrino spectrum as a Spectrum object.

        Args:
            cooling_time: The time in years since the cask was removed from the reactor.
                Note that this has to be greater than or equal to the
                initial_cooling_time of the cask.

        """
        spectra = self._get_component_spectra(cooling_time)

        # Equalise all the spectra to 1keV bins to allow combining,
        # going from 0 to the maximum energy across all spectra.
        max_energy = max(spec.energy[-1] for spec in spectra)
        for spec in spectra:
            spec.equalise(width=1, min_energy=0, max_energy=max_energy)

        # Sum all the spectra to get the total cask spectrum
        total_spec = spectra[0]
        for spec in spectra[1:]:
            total_spec = total_spec + spec
        total_spec.name = self.name
        return total_spec
