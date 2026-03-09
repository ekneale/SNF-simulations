"""Functions for loading and manipulating spectra."""

from pathlib import Path

import numpy as np

from .data import load_spectrum
from .utils import linear_interpolate_with_errors, sample_histogram


class Spectrum:
    """Class to represent an antineutrino spectrum.

    Attributes:
        energy: Array of energy values (keV), representing histogram bin edges.
        flux: Array of antineutrino flux values (keV^-1) for each energy bin.
            Note for N energy bin edges there should be N-1 flux values,
            with each pair of adjacent energy values defining the upper and lower
            edges of each bin.
        errors: Array of uncertainties for the flux values.
            Array length should be the same as flux.
        name: Name of the spectrum.

    """

    def __init__(
        self,
        energy: np.ndarray,
        flux: np.ndarray,
        errors: np.ndarray,
        name: str = "Spectrum",
    ) -> None:
        """Initialize the Spectrum object."""
        self.energy = energy
        self.flux = flux
        self.errors = errors
        self.name = name

        if self.energy.ndim != 1 or self.flux.ndim != 1 or self.errors.ndim != 1:
            msg = "Energy, flux and errors must be 1D arrays"
            raise ValueError(msg)
        if len(self.flux) != len(self.energy) - 1:
            msg = "Flux array must have length len(energy) - 1"
            raise ValueError(msg)
        if len(self.flux) != len(self.errors):
            msg = "Flux and errors arrays must have the same length"
            raise ValueError(msg)

    def __repr__(self) -> str:
        """Return a string representation of the Spectrum object."""
        try:
            repr_str = f'<Spectrum "{self.name}", '
            repr_str += f"energy_range=({self.energy[0]}-{self.energy[-1]} keV)>"
        except AttributeError:
            return "<Spectrum (uninitialized)>"
        else:
            return repr_str

    @classmethod
    def from_isotope(
        cls,
        name: str,
    ) -> "Spectrum":
        """Create a Spectrum object from an isotope name."""
        # The IAEA data files give equal arrays of energy, flux, and uncertainty.
        data = load_spectrum(name)
        energy_points, flux_points, error_points = data[:, 0], data[:, 1], data[:, 2]

        # Histogram representation requires N+1 edges for N bins.
        # The last energy point is used as the upper edge of the final bin,
        # so the last flux/error value is discarded.
        energy = energy_points
        flux = flux_points[:-1]
        errors = error_points[:-1]
        return cls(
            energy,
            flux,
            errors,
            name=name,
        )

    @classmethod
    def from_file(
        cls,
        filename: Path | str,
    ) -> "Spectrum":
        """Create a Spectrum object from a CSV file written using write_csv."""
        with open(filename) as f:
            header = f.readline().strip()
            # Name should have been saved in the header
            name = "Spectrum" if not header.startswith("#") else header[1:].strip()
        # The data should have columns: energy_lower, energy_upper, flux, error
        data = np.loadtxt(filename, delimiter=",", skiprows=2)
        lower_edges = data[:, 0]
        upper_edges = data[:, 1]
        energy = np.concatenate((lower_edges, [upper_edges[-1]]))
        flux = data[:, 2]
        errors = data[:, 3]
        return cls(
            energy,
            flux,
            errors,
            name=name,
        )

    def equalise(
        self,
        width: float = 1,
        min_energy: float | None = None,
        max_energy: float | None = None,
    ) -> None:
        """Convert the spectrum to have equal bin widths.

        Bins are spaced from min_energy to max_energy with the requested width.

        Args:
            width: Target bin width (keV). Must be positive.
            min_energy: Minimum energy for the new spectrum (keV).
                If None, uses the current minimum energy edge.
            max_energy: Maximum energy for the new spectrum (keV).
                If None, uses the current maximum energy edge.

        """
        if width <= 0:
            msg = "width must be a positive value"
            raise ValueError(msg)
        if min_energy is None:
            min_energy = float(self.energy[0])
        if max_energy is None:
            max_energy = float(self.energy[-1])
        if max_energy <= min_energy:
            msg = "max_energy must be greater than min_energy"
            raise ValueError(msg)
        if min_energy + width > max_energy:
            msg = "width is too large for the given energy range"
            raise ValueError(msg)

        # Interpolate to the new binning and propagate errors
        new_edges = np.arange(min_energy, max_energy + width, width)
        new_flux, new_errors = linear_interpolate_with_errors(
            self.energy,
            self.flux,
            self.errors,
            new_edges,
        )

        # Apply the new values to this Spectrum instance in place
        self.energy = new_edges
        self.flux = new_flux
        self.errors = new_errors

    def __add__(self, other: "Spectrum") -> "Spectrum":
        """Add another Spectrum to this one by summing the flux values.

        The energy bins of the two spectra must be the same.

        Args:
            other: Another Spectrum object to add to this one.

        Returns:
            A new Spectrum object representing the sum of the two spectra.

        """
        if not np.allclose(self.energy, other.energy):
            msg = "Energy bins of the two spectra must be the same to add them."
            raise ValueError(msg)
        new_flux = self.flux + other.flux
        new_errors = np.sqrt(self.errors**2 + other.errors**2)
        new_name = self.name + " + " + other.name
        return Spectrum(self.energy, new_flux, new_errors, name=new_name)

    def __mul__(self, factor: float) -> "Spectrum":
        """Multiply the spectrum by a scalar factor.

        Args:
            factor: The scaling factor to apply to the spectrum.

        Returns:
            A new Spectrum object representing the scaled spectrum.

        """
        new_flux = self.flux * factor
        new_errors = self.errors * abs(factor)
        return Spectrum(self.energy, new_flux, new_errors, name=self.name)

    def sample(self, samples: int = 100, seed: int | None = None) -> np.ndarray:
        """Sample the spectrum to simulate what a detector could observe.

        Args:
            samples: Number of samples to draw from the spectrum.
            seed: Random seed for reproducibility.


        Returns:
            Array of sampled energies.

        """
        return sample_histogram(self.energy, self.flux, samples, seed)

    def integrate(
        self,
        lower_energy: float | None = None,
        upper_energy: float | None = None,
    ) -> float:
        """Integrate the spectrum over an energy range.

        This provides the Spectrum-class equivalent of integrating a ROOT
        histogram for flux calculations.

        Args:
            lower_energy: Lower energy bound in keV.
                If None, uses the minimum energy from the spectrum.
            upper_energy: Upper energy bound in keV.
                If None, uses the maximum energy from the spectrum.

        Returns:
            Integrated spectrum value over the requested range.

        """
        if upper_energy is None:
            upper_energy = float(self.energy[-1])
        if lower_energy is None:
            lower_energy = float(self.energy[0])
        if upper_energy <= lower_energy:
            msg = "upper_energy must be greater than lower_energy"
            raise ValueError(msg)

        # We need to account for if the upper or lower bounds fall partially within a
        # bin instead of at a bin edge.
        # Calculate the overlap of each bin with the integration range,
        # and weight the flux in each bin by this overlap length.
        # For most bins this will be either 0 (outside the range) or the full bin width
        # (fully within the range), but for the two bins at the edges of the range
        # it could be a partial overlap.
        # Using np.clip gives an array where each value is the length of the overlap of
        # that bin with the integration range.
        # Then multiply the flux in each bin by the overlap length, and sum to get the
        # total integrated flux.
        lower_edges = self.energy[:-1]
        upper_edges = self.energy[1:]
        overlap_length = np.clip(
            np.minimum(upper_edges, upper_energy)
            - np.maximum(lower_edges, lower_energy),
            a_min=0,
            a_max=None,
        )
        return float(np.sum(self.flux * overlap_length))

    def write_csv(self, output_filename: Path | str = "") -> None:
        """Output energy and flux data to CSV file.

        Output format is:
            energy_lower, energy_upper, flux, error
            for each bin, where energy_lower and energy_upper are the lower and upper
            edges of the energy bin (in keV), flux is the flux value for that bin
            (in keV^-1), and error is the uncertainty for that flux value.

        Args:
            output_filename: The name of the output CSV file.

        """
        if not output_filename:
            output_filename = self.name.replace(" ", "_") + ".csv"
        if isinstance(output_filename, str) and not output_filename.endswith(".csv"):
            output_filename += ".csv"
        elif isinstance(output_filename, Path) and output_filename.suffix != ".csv":
            output_filename = output_filename.with_suffix(".csv")

        lower_edges = self.energy[:-1]
        upper_edges = self.energy[1:]
        data = np.column_stack((lower_edges, upper_edges, self.flux, self.errors))
        np.savetxt(
            output_filename,
            data,
            fmt=("%.1f", "%.1f", "%.6e", "%.6e"),
            delimiter=",",
            header=f"# {self.name}\nenergy_lower,energy_upper,flux,error",
            comments="",
        )
