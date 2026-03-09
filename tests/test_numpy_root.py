"""Temporary tests between legacy ROOT functions and the new Numpy-based classes.

These tests will be skipped if ROOT is not installed, and will be removed
once the ROOT-based code is fully removed.
"""

import numpy as np
import pytest

try:
    import ROOT
except ImportError:
    pytest.skip("skipping ROOT tests", allow_module_level=True)  # ty: ignore

from snf_simulations.cask import Cask
from snf_simulations.data import (
    load_antineutrino_data,
    load_isotope_data,
    load_reactor_data,
)
from snf_simulations.physics import DecayChain, get_decay_mass, get_isotope_activity
from snf_simulations.spec import Spectrum

from .test_spec import _mock_data

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts

ROOT.TH1.AddDirectory(False)  # Prevent ROOT from keeping histograms in memory


_MOLAR_MASSES, _ = load_isotope_data()
ISOTOPES = list(_MOLAR_MASSES.keys())
ISOTOPE_SPECTRA = load_antineutrino_data(ISOTOPES)


def _get_mock_spectra() -> tuple[Spectrum, ROOT.TH1D]:
    """Create a mock Spectrum and corresponding ROOT TH1D histogram."""
    # Use the same mock data as in test_spec.
    energy, flux, errors = _mock_data()

    # Create the Spectrum instance
    # Note we miss off the final flux and error values, but still pass the full
    # array of energy values.
    spec = Spectrum(energy, flux[:-1], errors[:-1], name="mock")

    # Create the ROOT histogram
    root_spec = _create_root_spec(energy, flux, errors, "mock")

    return spec, root_spec


def _root_to_arrays(spec: ROOT.TH1D) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert a ROOT TH1D histogram to Numpy arrays."""
    # ROOT bins are 1-indexed, so range from 1 to n_bins
    n_bins = spec.GetNbinsX()
    i_bins = range(1, n_bins + 1)

    # Get the edges for each bin from the ROOT histogram.
    # Note we have to specifically include the upper edge of the final bin.
    lower_edges = [spec.GetBinLowEdge(i) for i in i_bins]
    upper_edge = spec.GetBinLowEdge(n_bins) + spec.GetBinWidth(n_bins)
    edges = np.array([*lower_edges, upper_edge], dtype=float)

    # Get the content and errors from the ROOT histogram bins.
    content = np.array([spec.GetBinContent(i) for i in i_bins], dtype=float)
    errors = np.array([spec.GetBinError(i) for i in i_bins], dtype=float)

    # edges should be len(nbins + 1), content and errors should be len(n_bins)
    return edges, content, errors


def _compare_spectra(spec: Spectrum, root_spec: ROOT.TH1D) -> None:
    """Compare a Spectrum instance to a ROOT TH1D histogram."""
    root_energy, root_flux, root_errors = _root_to_arrays(root_spec)

    assert len(spec.energy) == len(root_energy), (
        f"Spectrum '{spec.name}' energy array has different length to ROOT energy array"
        f" ({len(spec.energy)} vs {len(root_energy)})"
    )
    assert len(spec.flux) == len(root_flux), (
        f"Spectrum '{spec.name}' flux array has different length to ROOT flux array"
        f" ({len(spec.flux)} vs {len(root_flux)})"
    )
    assert len(spec.errors) == len(root_errors), (
        f"Spectrum '{spec.name}' error array has different length to ROOT error array"
        f" ({len(spec.errors)} vs {len(root_errors)})"
    )

    assert np.allclose(spec.energy, root_energy), (
        f"Spectrum '{spec.name}' energy edges do not match ROOT edges"
    )
    assert np.allclose(spec.flux, root_flux), (
        f"Spectrum '{spec.name}' flux values do not match ROOT contents"
    )
    assert np.allclose(spec.errors, root_errors), (
        f"Spectrum '{spec.name}' error values do not match ROOT errors"
    )


def _create_root_spec(
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


def _equalise_root_spec(
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


def _scale_root_spec(
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


def _load_root_spec(  # noqa: PLR0913
    data: np.ndarray,
    name: str,
    mass: float,
    molar_mass: float,
    half_life: float,
    removal_time: float,
    max_energy: int | None = None,
    min_energy: int = 0,
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
    spec = _create_root_spec(data[:, 0], data[:, 1], data[:, 2], name)
    if max_energy is None:
        max_energy = int(np.ceil(max(data[:, 0])))
    spec_equal = _equalise_root_spec(spec, max_energy, min_energy)
    spec_scaled = _scale_root_spec(
        spec_equal, mass, molar_mass, half_life, removal_time
    )
    spec_scaled.SetTitle(name)
    return spec_scaled


def _add_root_spec(spectra: ROOT.TList) -> ROOT.TH1D:
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


def _get_total_root_spec(
    cask_name: str,
    isotope_proportions: dict,
    total_mass: float = 1000,
    removal_time: float = 0,
    max_energy: int | None = None,
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
        spec = _load_root_spec(
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
            spec = _load_root_spec(
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
    total_spec = _add_root_spec(spectra)
    total_spec.SetTitle("Total Spectrum")
    if max_energy is not None:
        total_spec.GetXaxis().SetRangeUser(0, max_energy)
    return total_spec


def test_spectrum() -> None:
    """Test that the Spectrum class matches ROOT."""
    spec, root_spec = _get_mock_spectra()

    # Compare the arrays to ensure they match
    _compare_spectra(spec, root_spec)


def test_equalise() -> None:
    """Test that Spectrum.equalise matches ROOT."""
    spec, root_spec = _get_mock_spectra()

    # Equalise the spectra
    spec.equalise(width=1, min_energy=0, max_energy=3)
    root_equal = _equalise_root_spec(root_spec, max_energy=3, min_energy=0)

    # Compare the arrays to ensure they match
    _compare_spectra(spec, root_equal)


def test_scale() -> None:
    """Test that Spectrum scaling matches ROOT."""
    spec, root_spec = _get_mock_spectra()

    # Scale the spectra
    mass = 1000
    molar_mass = 100
    half_life = 30
    removal_time = 5
    activity = get_isotope_activity(mass, molar_mass, half_life, removal_time)
    spec_scaled = spec * activity
    root_scaled = _scale_root_spec(root_spec, mass, molar_mass, half_life, removal_time)

    # Compare the arrays to ensure they match
    _compare_spectra(spec_scaled, root_scaled)


def test_add() -> None:
    """Test that adding two Spectra matches adding two ROOT spectra."""
    spec1, root_spec1 = _get_mock_spectra()
    spec2, root_spec2 = _get_mock_spectra()

    # Add the spectra together
    spec_sum = spec1 + spec2
    root_spectra = ROOT.TList()
    root_spectra.Add(root_spec1)
    root_spectra.Add(root_spec2)
    root_sum = _add_root_spec(root_spectra)

    # Compare the arrays to ensure they match
    _compare_spectra(spec_sum, root_sum)


def test_sample() -> None:
    """Test that Spectrum.sample produces a distribution comparable to ROOT's."""
    spec, root_spec = _get_mock_spectra()

    # Sample from both spectra
    samples = 100000
    spec_samples = spec.sample(samples=samples)
    root_samples = np.array([root_spec.GetRandom() for _ in range(samples)])
    assert len(spec_samples) == samples, (
        "Spectrum.sample did not produce the expected number of samples"
    )
    assert len(root_samples) == samples, (
        "ROOT.GetRandom did not produce the expected number of samples"
    )

    # We can't compare the samples directly due to the random selection, but we can
    # compare their distributions to the expected probabilities from the histogram.
    # First, compare against the expected probabilities which are proportional to
    # the (normalised) flux values in each bin.
    expected_probabilities = spec.flux / np.sum(spec.flux)
    spec_counts, _ = np.histogram(spec_samples, bins=spec.energy)
    root_counts, _ = np.histogram(root_samples, bins=spec.energy)

    # We allow a reasonable 5-sigma error margin, where sigma is the standard deviation
    # of the binomial distribution for each bin.
    sigma = np.sqrt(expected_probabilities * (1 - expected_probabilities) / samples)
    spec_probabilities = spec_counts / samples
    root_probabilities = root_counts / samples
    assert np.all(np.abs(spec_probabilities - expected_probabilities) < 5 * sigma), (
        "Spectrum.sample bin probabilities deviate more than expected from "
        "the target histogram distribution"
    )
    assert np.all(np.abs(root_probabilities - expected_probabilities) < 5 * sigma), (
        "ROOT.GetRandom bin probabilities deviate more than expected from "
        "the target histogram distribution"
    )

    # Two-sample chi-squared consistency check between the two outputs.
    combined_counts = root_counts + spec_counts
    # Exclude bins where both counts are zero to avoid division by zero in chi2.
    valid = combined_counts > 0
    chi2 = np.sum(
        (root_counts[valid] - spec_counts[valid]) ** 2 / combined_counts[valid],
    )
    degrees_of_freedom = int(np.count_nonzero(valid) - 1)
    assert degrees_of_freedom > 0

    # The standard deviation of the chi-squared distribution is sqrt(2*ndf),
    # again allow a 5-sigma deviation from the mean.
    sigma = np.sqrt(2 * degrees_of_freedom)
    assert chi2 < degrees_of_freedom + 5 * sigma, (
        f"ROOT and Spectrum sampled distributions differ too much: "
        f"chi2={chi2:.3f}, ndf={degrees_of_freedom}"
    )


def test_integrate() -> None:
    """Test that Spectrum.integrate works for arbitrary (non-equalised) bins."""
    spec, root_spec = _get_mock_spectra()

    # Test integration over full range.
    # Flux is dN/dE (keV^-1), so integral is sum(flux * bin_width).
    expected_integral = np.sum(spec.flux * np.diff(spec.energy))
    spec_integral = spec.integrate(spec.energy[0], spec.energy[-1])
    root_integral = root_spec.Integral(1, root_spec.GetNbinsX(), "width")

    assert np.isclose(spec_integral, expected_integral), (
        f"Spectrum.integrate({spec.energy[0]}, {spec.energy[-1]}) = {spec_integral}"
        f" does not match expected value {expected_integral}"
    )
    assert np.isclose(root_integral, expected_integral), (
        f"ROOT spec.Integral(1, {root_spec.GetNbinsX()}, 'width') = {root_integral}"
        f" does not match expected value {expected_integral}"
    )

    # Test integration over partial range (still full bins).
    lower_energy = 0.5
    upper_energy = 3.0
    lower_index = 1
    upper_index = 4
    expected_integral = np.sum(
        spec.flux[lower_index:upper_index]
        * np.diff(spec.energy)[lower_index:upper_index],
    )
    spec_integral = spec.integrate(lower_energy, upper_energy)
    root_integral = root_spec.Integral(lower_index + 1, upper_index, "width")
    assert np.isclose(spec_integral, expected_integral), (
        f"Spectrum.integrate({lower_energy}, {upper_energy}) = {spec_integral}"
        f" does not match expected value {expected_integral}"
    )
    assert np.isclose(root_integral, expected_integral), (
        f"ROOT spec.Integral({lower_index + 1}, {upper_index}, 'width') = "
        f"{root_integral}"
        f" does not match expected value {expected_integral}"
    )

    # Currently, we only ever use ROOT's Integral on spectra that have been equalised
    # to 1 keV bins from 0 to 6000. That means we don't have to worry about partial
    # bin widths and can use the energy range values as bin values.
    # For the Spectrum class we want to be able to integrate over arbitrary energy
    # ranges, including partial bins.
    # We can test that elsewhere, since trying to get ROOT to do partial bin integration
    # is non-trivial.


@pytest.mark.parametrize("isotope", ISOTOPES)
def test_spec_real(isotope: str) -> None:
    """Test that Spectra from real data matches ROOT."""
    # Create the Spectrum instance and equalise it
    spec = Spectrum.from_isotope(isotope)
    spec.equalise(width=1, min_energy=0)

    # Create the ROOT spectrum and equalise it
    isotope_data = ISOTOPE_SPECTRA[isotope]
    energy = isotope_data[:, 0]
    flux = isotope_data[:, 1]
    errors = isotope_data[:, 2]
    root_spec = _create_root_spec(energy, flux, errors, isotope)
    max_energy = int(np.ceil(max(isotope_data[:, 0])))
    root_equal = _equalise_root_spec(root_spec, max_energy=max_energy, min_energy=0)

    # Compare the arrays to ensure they match
    _compare_spectra(spec, root_equal)


@pytest.mark.parametrize("total_mass,removal_time", [(1000, 0), (5000, 10)])
def test_cask(total_mass: float, removal_time: float) -> None:
    """Test that the Cask class matches ROOT."""
    # Test with only two isotopes
    # Note for removal_time>0 Sr90 will decay into Y90
    isotope_proportions = {"Sr90": 0.5, "Cs137": 0.5}

    # Create the Cask and get the total spectra at the given removal time
    cask = Cask(
        isotope_proportions=isotope_proportions,
        total_mass=total_mass,
    )
    spec = cask.get_total_spectrum(removal_time=removal_time)

    # Do the same with the ROOT function
    root_spec = _get_total_root_spec(
        cask_name="Test",
        isotope_proportions=isotope_proportions,
        total_mass=total_mass,
        removal_time=removal_time,
    )

    # Compare the spectra to ensure they match
    _compare_spectra(spec, root_spec)


@pytest.mark.parametrize(
    "reactor,total_mass,removal_time",
    [
        ("sizewell", 1000, 0),
        ("sizewell", 5000, 10),
        ("hartlepool", 1000, 0),
        ("hartlepool", 5000, 10),
    ],
)
def test_cask_real(reactor: str, total_mass: float, removal_time: float) -> None:
    """Test that the Cask class from real data matches ROOT."""
    isotope_proportions = load_reactor_data(reactor)

    # Create the Cask and get the total spectra at the given removal time
    cask = Cask(
        isotope_proportions=isotope_proportions,
        total_mass=total_mass,
    )
    spec = cask.get_total_spectrum(removal_time=removal_time)

    # Do the same with the ROOT function
    root_spec = _get_total_root_spec(
        cask_name="Test",
        isotope_proportions=isotope_proportions,
        total_mass=total_mass,
        removal_time=removal_time,
    )

    # Compare the spectra to ensure they match
    _compare_spectra(spec, root_spec)
