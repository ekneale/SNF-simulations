"""Functions for loading and manipulating spectra."""

import numpy as np
import ROOT

from .physics import get_isotope_activity


def create_spec(energy, dn_de, errors, name):
    """Load spectrum data into a ROOT TH1D histogram.

    Args:
        energy (array-like): Array of energy bin edges (keV).
        dn_de (array-like): Array of dN/dE values corresponding to each energy bin.
        errors (array-like): Array of errors for each dN/dE value.
        name (str): Name of the histogram.

    Returns:
        ROOT.TH1D: Histogram representing the spectrum.

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
    for e, dn in zip(energy, dn_de):
        spec.Fill(e, dn)
    for i, err in enumerate(errors, start=1):
        spec.SetBinError(i, err)

    # Formatting for display
    spec.SetStats(0)
    spec.GetXaxis().SetTitle("Energy [keV]")
    spec.GetYaxis().SetTitle("#frac{dN}{dE} [keV^{-1}]")

    return spec


def equalise_spec(spec, max_energy, min_energy=0):
    """Convert a ROOT TH1D histogram to have equal bin widths.

    Each bin will be 1 keV wide, spaced from min_energy to max_energy.

    Args:
        spec (ROOT.TH1D): Original histogram with variable bin widths.
        max_energy (float): Maximum energy for the new histogram.
        min_energy (float, default=0): Minimum energy for the new histogram.

    Returns:
        ROOT.TH1D: New histogram with equal bin widths.

    """
    # Create new bin edges and new bin centres.
    new_edges = np.linspace(min_energy, max_energy, (max_energy - min_energy) + 1)
    new_centres = new_edges[:-1] + 0.5

    # Interpolate the content of the original histogram to the new bin centres
    new_content = [spec.Interpolate(centre) for centre in new_centres]

    # Calculate the new errors for interpolated bins
    new_errors = []
    for centre in new_centres:
        # The TH1D::Interpolate function looks for the two closest bins surrounding the
        # given point.
        # If it's below the first bin centre or above the final bin centre it just
        # returns the content of the first/last bin, so here we do the same.
        if centre < spec.GetBinCenter(1):
            new_errors.append(spec.GetBinError(1))
            continue
        if centre >= spec.GetBinCenter(spec.GetNbinsX()):
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
    for centre, content in zip(new_centres, new_content):
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


def scale_spec(spec, mass, molar_mass, half_life, removal_time):
    """Scale a ROOT TH1D histogram spectrum by the activity of the isotope.

    Args:
        spec (ROOT.TH1D): Histogram representing the spectrum to be scaled.
        mass (float): Mass of the isotope in kg.
        molar_mass (float): Molar mass of the isotope in g/mol.
        half_life (float): Half-life of the isotope in years.
        removal_time (float): Time since removal from reactor in years.

    Returns:
        ROOT.TH1D: Scaled histogram.

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


def load_spec(
    data,
    name,
    mass,
    molar_mass,
    half_life,
    removal_time,
    max_energy=None,
    min_energy=0,
):
    """Load, equalise and scale a spectrum from data.

    Combines the create_spec, equalise_spec and scale_spec functions.
    """
    spec = create_spec(data[:, 0], data[:, 1], data[:, 2], name)
    if max_energy is None:
        max_energy = int(np.floor(max(data[:, 0])))
    spec_equal = equalise_spec(spec, max_energy, min_energy)
    spec_scaled = scale_spec(spec_equal, mass, molar_mass, half_life, removal_time)
    spec_scaled.SetTitle(name)
    return spec_scaled


def add_spec(spectra):
    """Combine a ROOT.TList of spectra into one total spectrum."""
    total_spec = spectra[0].Clone("combined")
    total_spec.Reset()
    total_spec.Merge(spectra)
    return total_spec


def write_spec(spec, output_filename):
    """Output energy and flux data to CSV file.

    Args:
        spec (ROOT.TH1D): The flux spectrum histogram.
        output_filename (str): The name of the output CSV file.

    Returns:
        data: A 2D numpy array with energy (keV) and flux (keV^{-1} s^{-1}) columns.

    """
    # Extract energy and flux data from the histogram
    n_bins = spec.GetNbinsX()
    energy = np.array([spec.GetBinCenter(i) for i in range(1, n_bins + 1)])
    flux = np.array([spec.GetBinContent(i) for i in range(1, n_bins + 1)])
    data = np.column_stack((energy, flux))

    # Save energy and flux data as a csv file
    if not output_filename.endswith(".csv"):
        output_filename += ".csv"
    np.savetxt(
        output_filename,
        data,
        fmt=("%.1f", "%.6e"),
        delimiter=",",
        header="energy,flux",
        comments="",
    )

    return data
