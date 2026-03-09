"""Functions for loading and manipulating spectra."""

from pathlib import Path

import numpy as np
import ROOT

from .data import load_spectrum
from .physics import get_isotope_activity
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


def create_spec(
    energy: np.ndarray,
    flux: np.ndarray,
    errors: np.ndarray,
    name: str,
) -> ROOT.TH1D:
    """Load spectrum data into a ROOT TH1D histogram.

    Args:
        energy: Array of energy bin edges (keV).
        flux: Array of flux values (dN/dE) corresponding to each energy bin (keV^-1).
        errors: Array of errors for each flux value.
        name: Name of the histogram.

    Returns:
        Histogram representing the spectrum.

    """
    # Create a histogram, using the energies as the bin edges.
    # Arguments are:
    #   name
    #   histogram title
    #   number of bins
    #   array of lower bin edges
    # Note that number of bins is len(energy) - 1 as there is one less bin than edges,
    # and the last energy value is assumed to be the upper edge of the final bin.
    spec = ROOT.TH1D(name, name, len(energy) - 1, np.array(energy))

    # Fill the bins
    # Note that ROOT histograms are 1-indexed, so bin 1 takes array value [0] etc.
    # For spec.Fill this doesn't matter, as it fills each bin sequentially.
    # But for spec.SetBinError we loop starting at bin 1, so need the error from [i-1].
    # Using enumerate() with start=1 will handle this nicely.
    for e, f in zip(energy, flux, strict=True):
        spec.Fill(e, f)
    for i, err in enumerate(errors, start=1):
        spec.SetBinError(i, err)

    # Formatting for display
    spec.SetStats(0)
    spec.GetXaxis().SetTitle("Energy [keV]")
    spec.GetYaxis().SetTitle("#frac{dN}{dE} [keV^{-1}]")

    return spec


def equalise_spec(
    spec: ROOT.TH1D,
    max_energy: int,
    min_energy: int = 0,
) -> ROOT.TH1D:
    """Convert a ROOT TH1D histogram to have equal bin widths.

    Each bin will be 1 keV wide, spaced from min_energy to max_energy.

    Args:
        spec: Original histogram with variable bin widths.
        max_energy: Maximum energy for the new histogram.
        min_energy: Minimum energy for the new histogram.

    Returns:
        New histogram with equal bin widths.

    """
    # Create new bin edges and new bin centres.
    new_edges = np.linspace(min_energy, max_energy, (max_energy - min_energy) + 1)
    new_centres = new_edges[:-1] + 0.5

    first_edge = spec.GetBinLowEdge(1)
    n_bins = spec.GetNbinsX()
    last_edge = spec.GetBinLowEdge(n_bins) + spec.GetBinWidth(n_bins)
    lower_edges = new_edges[:-1]
    upper_edges = new_edges[1:]
    outside = (upper_edges <= first_edge) | (lower_edges >= last_edge)

    # Interpolate the content of the original histogram to the new bin centres.
    # Outside the source range, bins are set to zero.
    new_content = [
        0.0 if is_outside else spec.Interpolate(centre)
        for centre, is_outside in zip(new_centres, outside, strict=True)
    ]

    # Calculate the new errors for interpolated bins
    new_errors = []
    for centre, is_outside in zip(new_centres, outside, strict=True):
        # The TH1D::Interpolate function looks for the two closest bins surrounding the
        # given point.
        # Outside the original range, extrapolated bins are set to zero.
        if is_outside:
            new_errors.append(0.0)
            continue
        if centre <= spec.GetBinCenter(1) or np.isclose(centre, spec.GetBinCenter(1)):
            new_errors.append(spec.GetBinError(1))
            continue
        if centre >= spec.GetBinCenter(spec.GetNbinsX()) or np.isclose(
            centre, spec.GetBinCenter(spec.GetNbinsX())
        ):
            new_errors.append(spec.GetBinError(spec.GetNbinsX()))
            continue

        # Find which of the old bins the new centre would fall within.
        idx = spec.FindBin(centre)

        # Find which of the surrounding bins is lower and which is upper.
        # We need to check which side of the centre of the found bin the new centre is.
        if centre < spec.GetBinCenter(idx):
            lower_centre = spec.GetBinCenter(idx - 1)
            upper_centre = spec.GetBinCenter(idx)
            lower_error = spec.GetBinError(idx - 1)
            upper_error = spec.GetBinError(idx)
        else:
            lower_centre = spec.GetBinCenter(idx)
            upper_centre = spec.GetBinCenter(idx + 1)
            lower_error = spec.GetBinError(idx)
            upper_error = spec.GetBinError(idx + 1)

        # Linear interpolation of the new error.
        # We find the new error by scaling the errors of the surrounding two bins
        # by the distance of the new centre to each of them.
        dist = abs(upper_centre - lower_centre)
        dist_lower = abs(centre - lower_centre)
        dist_upper = abs(upper_centre - centre)
        weight_lower = (dist - dist_lower) / dist
        weight_upper = (dist - dist_upper) / dist
        error = np.sqrt(
            (weight_lower * lower_error) ** 2 + (weight_upper * upper_error) ** 2,
        )
        new_errors.append(error)

    # Create the new histogram
    name = spec.GetName()
    title = spec.GetTitle()
    spec_equal = ROOT.TH1D(
        str(name),
        str(title) + " equal bin widths",
        len(new_centres),
        new_edges,
    )

    # Fill the new histogram with the interpolated content
    for centre, content in zip(new_centres, new_content, strict=True):
        spec_equal.Fill(centre, content)
    for i, err in enumerate(new_errors, start=1):
        spec_equal.SetBinError(i, err)

    # Formatting for display
    spec_equal.SetStats(0)
    spec_equal.GetXaxis().SetTitle("Energy [keV]")
    spec_equal.GetYaxis().SetTitle("#frac{dN}{dE} [keV^{-1}]")
    spec_equal.GetXaxis().SetLabelSize(0.05)
    spec_equal.GetYaxis().SetLabelSize(0.05)
    spec_equal.GetXaxis().SetTitleSize(0.047)
    spec_equal.GetYaxis().SetTitleSize(0.047)

    return spec_equal


def scale_spec(
    spec: ROOT.TH1D,
    mass: float,
    molar_mass: float,
    half_life: float,
    removal_time: float,
) -> ROOT.TH1D:
    """Scale a ROOT TH1D histogram spectrum by the activity of the isotope.

    Args:
        spec: Histogram representing the spectrum to be scaled.
        mass: Mass of the isotope in kg.
        molar_mass: Molar mass of the isotope in g/mol.
        half_life: Half-life of the isotope in years.
        removal_time: Time since removal from reactor in years.

    Returns:
        Scaled histogram.

    """
    # Calculate the activity of the isotope
    activity = get_isotope_activity(
        mass,
        molar_mass,
        half_life,
        removal_time,
    )

    # Scale the spectrum
    spec.Scale(activity)

    # Formatting for display
    spec.SetStats(0)
    title = spec.GetTitle()
    spec.SetTitle(title + " scaled")
    spec.GetXaxis().SetTitle("Energy [keV]")
    spec.GetYaxis().SetTitle("Relative Flux [keV^{-1}s^{-1}]")

    return spec


def load_spec(  # noqa: PLR0913
    data: np.ndarray,
    name: str,
    mass: float,
    molar_mass: float,
    half_life: float,
    removal_time: float,
    max_energy: float | None = None,
    min_energy: float = 0,
) -> ROOT.TH1D:
    """Load, equalise and scale a spectrum from data.

    Combines the create_spec, equalise_spec and scale_spec functions to load
    raw spectrum data, convert it to equal bin widths, and scale by isotope activity.

    Args:
        data: Array of spectrum data with columns [energy, flux, uncertainty].
        name: Name of the histogram.
        mass: Mass of the isotope in kg.
        molar_mass: Molar mass of the isotope in g/mol.
        half_life: Half-life of the isotope in years.
        removal_time: Time since removal from reactor in years.
        max_energy: Maximum energy for the new histogram (keV).
            If None, uses the maximum energy from the input data.
        min_energy: Minimum energy for the new histogram (keV).

    Returns:
        Loaded, equalised and scaled spectrum histogram.

    """
    # TODO: there are too many arguments to this function, but if ROOT is removed
    # and atomic properties can be fetched through mendeleev then it can
    # be simplified a lot.
    spec = create_spec(data[:, 0], data[:, 1], data[:, 2], name)
    if max_energy is None:
        max_energy = int(np.ceil(max(data[:, 0])))
    spec_equal = equalise_spec(spec, max_energy, min_energy)
    spec_scaled = scale_spec(spec_equal, mass, molar_mass, half_life, removal_time)
    spec_scaled.SetTitle(name)
    return spec_scaled


def add_spec(spectra: ROOT.TList) -> ROOT.TH1D:
    """Combine a ROOT.TList of spectra into one total spectrum.

    Merges multiple ROOT histograms into a single combined spectrum by summing
    the bin contents across all input histograms.

    Args:
        spectra: A ROOT TList containing ROOT.TH1D spectrum histograms.

    Returns:
        Combined spectrum with summed bin contents from all input spectra.

    """
    total_spec = spectra[0].Clone("combined")
    total_spec.Reset()
    total_spec.Merge(spectra)
    return total_spec


def sample_spec(
    spec: ROOT.TH1D,
    samples: int = 100,
    output_filename: Path | str | None = None,
) -> list[float]:
    """Sample the total spectrum to simulate what a detector could observe.

    Samples can be saved as CSV for further analysis.

    Args:
        spec: The flux spectrum histogram.
        samples: Number of samples to draw from the spectrum.
        output_filename: Filename to save the sampled data as a CSV file.
            If None, data is not saved.

    Returns:
        List of sampled flux values.

    """
    sampled_flux = [spec.GetRandom() for _ in range(samples)]

    if output_filename:
        # Save samples as a CSV file
        if isinstance(output_filename, str) and not output_filename.endswith(".csv"):
            output_filename += ".csv"
        elif isinstance(output_filename, Path) and output_filename.suffix != ".csv":
            output_filename = output_filename.with_suffix(".csv")
        np.savetxt(
            output_filename,
            sampled_flux,
            delimiter=",",
            comments="",
        )

    return sampled_flux


def write_spec(
    spec: ROOT.TH1D,
    output_filename: Path | str,
) -> np.ndarray:
    """Output energy and flux data to CSV file.

    Args:
        spec: The flux spectrum histogram.
        output_filename: The name of the output CSV file.

    Returns:
        A 2D numpy array with energy (keV) and flux (keV^{-1} s^{-1}) columns.

    """
    # Extract energy and flux data from the histogram
    n_bins = spec.GetNbinsX()
    energy = np.array([spec.GetBinCenter(i) for i in range(1, n_bins + 1)])
    flux = np.array([spec.GetBinContent(i) for i in range(1, n_bins + 1)])
    data = np.column_stack((energy, flux))

    # Save energy and flux data as a csv file
    if isinstance(output_filename, str) and not output_filename.endswith(".csv"):
        output_filename += ".csv"
    elif isinstance(output_filename, Path) and output_filename.suffix != ".csv":
        output_filename = output_filename.with_suffix(".csv")
    np.savetxt(
        output_filename,
        data,
        fmt=("%.1f", "%.6e"),
        delimiter=",",
        header="energy,flux",
        comments="",
    )

    return data
