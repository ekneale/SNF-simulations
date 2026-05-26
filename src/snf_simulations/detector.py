"""Functions for simulating antineutrino detectors."""

from .physics import calculate_flux_at_distance
from .spec import Spectrum


class Detector:
    """Class representing an antineutrino detector.

    Attributes:
        volume: Volume of the detector in cubic meters.
        proton_density: Proton density of the detector material per cubic centimeter.
        name: An optional name for the detector.

    """

    def __init__(
        self,
        volume: float,
        proton_density: float,
        name: str | None = None,
    ) -> None:
        """Initialize the Detector object."""
        self.volume = volume
        self.proton_density = proton_density
        self.number_of_protons = (self.volume * 1e6) * self.proton_density
        self.name = name

        if self.volume <= 0:
            msg = "Detector volume must be a positive value"
            raise ValueError(msg)
        if self.proton_density <= 0:
            msg = "Proton density must be a positive value"
            raise ValueError(msg)

    def __repr__(self) -> str:
        """Return a string representation of the Detector object."""
        try:
            name_str = f' "{self.name}"' if self.name is not None else ""
            repr_str = (
                f"<Detector{name_str}: "
                f"volume={self.volume:.3e} m3, "
                f"proton_density={self.proton_density:.3e} protons/cm3>"
            )
        except AttributeError:
            return "<Detector (uninitialized)>"
        else:
            return repr_str

    def calculate_event_rate(
        self,
        spec: Spectrum,
        distance: float,
        efficiency: float = 1,
    ) -> float:
        """Calculate the expected event rate in the detector for a given spectrum.

        Args:
            spec: Antineutrino spectrum.
            distance: Distance from the source to the detector in meters.
            efficiency: Detection efficiency.

        Returns:
            Event rate in s^-1.

        """
        # The total flux is the integral of the spectrum over an energy range.
        # For inverse beta decay the threshold is 1.806 MeV, so we integrate above this
        # energy to get the total flux of antineutrinos that can be detected.
        threshold_energy = 1806  # KeV
        if spec.energy[-1] < threshold_energy:
            msg = "Spectrum energy range does not cover the IBD threshold of 1.806 MeV"
            raise ValueError(msg)
        total_flux = spec.integrate(lower_energy=threshold_energy)

        # Calculate the antineutrino flux at a given distance from the source.
        flux_at_distance = calculate_flux_at_distance(total_flux, distance)

        # Calculate the event rate using the flux and number of target protons.
        cross_section = 1e-44  # cm^2, approximate IBD cross-section
        event_rate = self.number_of_protons * flux_at_distance * cross_section

        # Apply detection efficiency and return
        return event_rate * efficiency
