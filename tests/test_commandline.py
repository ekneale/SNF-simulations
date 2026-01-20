"""Basic functional tests to check the output of examples from the command line script.

These aren't 'real' tests yet, there area few issues like relative paths to data files
that need to be resolved before they can be fully integrated into a test suite.

For now this script has to be run from within the snf_simulations directory.

But at least it should check that the output from the command line script hasn't changed,
and give a baseline for future tests.
"""

import os

import numpy as np
import ROOT
from snf_simulations import flux, plotting

ROOT.TH1.AddDirectory(False)  # Prevent ROOT from keeping histograms in memory


def _load_output(filename):
    """Load data from an output file.

    The 'CSV' files written by the `flix.write_spec_single` and
    `flix.write_spec_multiple` functions are not really proper CSV files,
    so we have to do a bit of manual parsing here to get the data out.
    """
    with open(filename) as f:
        energy_line = f.readline()
        flux_line = f.readline()

    # Line format is:
    # "energy": [0.0005, 0.0015, ... 5.2995, 5.3005],\n
    # "(flux)": [58875835323384.055, ... 0.0, 0.0],\n
    energy = [float(e) for e in energy_line[11:-3].split(", ")]
    flux = [float(fx) for fx in flux_line[11:-3].split(", ")]

    return energy, flux


def _test_single_cask(site="sizewell", removal_times=[0.5, 1, 5, 10, 20]):
    """Test single cask spectrum output."""
    print(f"--- Testing single cask for {site.capitalize()} ---")

    # Create the spectrum plot
    # This function will pop up a ROOT window, then save the plot to a pdf.
    # What we really need though is the returned spectrum data, which will
    # be a single TH1D object. It's currently hardcoded to return the first from
    # the list of removal times, i.e. 0.5 years by default.
    spec_single = plotting.plot_single_cask(
        removal_times,
        Sizewell=site == "sizewell",
        HartlePool=site == "hartlepool",
    )

    # We can't easily check the output plot, but we can check the basics of the returned spectrum.
    assert isinstance(spec_single, ROOT.TH1D), "Returned spectrum is not a TH1D"
    counts, values = spec_single.counts(), spec_single.values()
    assert isinstance(counts, np.ndarray), "Counts are not a numpy array"
    assert isinstance(values, np.ndarray), "Values are not a numpy array"
    assert len(counts) == len(values) == 5312, "Data has wrong length"

    # Write the spectrum to a file and get back the energy and flux arrays
    energy_single, flux_single = flux.write_spec_single(
        spec_single,
        Sizewell=site == "sizewell",
        Hartlepool=site == "hartlepool",
    )

    # Check the returned arrays
    assert isinstance(energy_single, list), "Energy is not a list"
    assert isinstance(flux_single, list), "Flux is not a list"
    assert len(energy_single) == len(flux_single) == 5312, "Data has wrong length"

    # Load the data back from the output file
    # (this assumes the file has been written in the current working directory)
    filename = f"{site.capitalize()}_single_0.5.csv"
    energy_loaded, flux_loaded = _load_output(filename)
    assert energy_loaded == energy_single, "Loaded energy does not match"
    assert flux_loaded == flux_single, "Loaded flux does not match"

    # Compare to the reference data file
    # (this again assumes we're in the main package directory)
    ref_filename = f"../../tests/test_data/{site.capitalize()}_single_0.5.csv"
    energy_ref, flux_ref = _load_output(ref_filename)
    assert energy_ref == energy_single, "Reference energy does not match"
    assert flux_ref == flux_single, "Reference flux does not match"

    # Delete the output files
    os.remove(filename)
    os.remove(f"{site.capitalize()}_Spectra_0.5.pdf")

    print("--- Single tests passed ---")


def _test_multiple_casks(site="sizewell"):
    """Test multiple cask spectrum output."""
    print(f"--- Testing multiple casks for {site.capitalize()} ---")

    # Create the spectrum plot
    # This function actually doesn't create a plot, it just returns the
    # combined spectrum.
    # Also both functions have the removal times hardcoded, so we pass None for now.
    if site == "sizewell":
        spec_multiple = plotting.plot_multiple_casks_sizewell(removal_times=None)
    elif site == "hartlepool":
        spec_multiple = plotting.plot_multiple_casks_hartlepool(removal_times=None)

    # Check the returned spectrum.
    assert isinstance(spec_multiple, ROOT.TH1D), "Returned spectrum is not a TH1D"
    counts, values = spec_multiple.counts(), spec_multiple.values()
    assert isinstance(counts, np.ndarray), "Counts are not a numpy array"
    assert isinstance(values, np.ndarray), "Values are not a numpy array"
    assert len(counts) == len(values) == 5312, "Data has wrong length"

    # Write the spectrum to a file and get back the energy and flux arrays
    energy_multiple, flux_multiple = flux.write_spec_multiple(
        spec_multiple,
        Sizewell=site == "sizewell",
        Hartlepool=site == "hartlepool",
    )

    # Check the returned arrays
    assert isinstance(energy_multiple, list), "Energy is not a list"
    assert isinstance(flux_multiple, list), "Flux is not a list"
    assert len(energy_multiple) == len(flux_multiple) == 5312, "Data has wrong length"

    # Load the data back from the output file
    # (this assumes the file has been written in the current working directory)
    filename = f"{site.capitalize()}_multiple.csv"
    energy_loaded, flux_loaded = _load_output(filename)
    assert energy_loaded == energy_multiple, "Loaded energy does not match"
    assert flux_loaded == flux_multiple, "Loaded flux does not match"

    # Compare to the reference data file
    # (this again assumes we're in the main package directory)
    ref_filename = f"../../tests/test_data/{site.capitalize()}_multiple.csv"
    energy_ref, flux_ref = _load_output(ref_filename)
    assert energy_ref == energy_multiple, "Reference energy does not match"
    assert flux_ref == flux_multiple, "Reference flux does not match"

    # Delete the output file
    os.remove(filename)

    print("--- Multiple tests passed ---")


def _test_flux_calc(site="sizewell", removal_times=[0.5, 1, 5, 10, 20]):
    """Test flux calculation at a given distance."""
    print(f"--- Testing flux calculation for {site.capitalize()} ---")

    # Get the single and multiple cask spectra
    # Annoyingly this will pop up the ROOT window and save the pdf again.
    spec_single = plotting.plot_single_cask(
        removal_times,
        Sizewell=site == "sizewell",
        HartlePool=site == "hartlepool",
    )
    if site == "sizewell":
        spec_multiple = plotting.plot_multiple_casks_sizewell(removal_times=None)
    elif site == "hartlepool":
        spec_multiple = plotting.plot_multiple_casks_hartlepool(removal_times=None)
    os.remove(f"{site.capitalize()}_Spectra_0.5.pdf")

    # Calculate the single cask flux at 40m
    # Returns a single value of flux cm^-2 s^-1
    # The function will also print out different conversions and rate limits,
    # but they don't get returned.
    flux_single_40 = flux.flux_calc(spec_single, distance_m=40)
    assert isinstance(flux_single_40, float), "Single flux is not a float"
    if site == "sizewell":
        assert flux_single_40 == 11992567783.00658, "Single flux value does not match"
    else:
        assert flux_single_40 == 4942443091.633026, "Single flux value does not match"

    # Calculate the multiple cask flux at 40m
    flux_multiple_40 = flux.flux_calc(spec_multiple, distance_m=40)
    assert isinstance(flux_multiple_40, float), "Multiple flux is not a float"
    if site == "sizewell":
        assert flux_multiple_40 == 13067897209.144945, (
            "Multiple flux value does not match"
        )
    else:
        assert flux_multiple_40 == 1443268305.2325196, (
            "Multiple flux value does not match"
        )

    print("--- Flux calculation tests passed ---")


def _test_multiple_plot(site="sizewell", removal_times=[0.5, 1, 5, 10, 20]):
    """Test plotting multiple flux spectra on one graph."""
    print(f"--- Testing multiple flux plot for {site.capitalize()} ---")

    # Get the single and multiple cask spectra
    # Annoyingly this will pop up the ROOT window and save the pdf again.
    spec_single = plotting.plot_single_cask(
        removal_times,
        Sizewell=site == "sizewell",
        HartlePool=site == "hartlepool",
    )
    if site == "sizewell":
        spec_multiple = plotting.plot_multiple_casks_sizewell(removal_times=None)
        os.remove("Sizewell_Spectra_0.5.pdf")
    elif site == "hartlepool":
        spec_multiple = plotting.plot_multiple_casks_hartlepool(removal_times=None)
        os.remove("Hartlepool_Spectra_0.5.pdf")

    # Write out the single and multiple spectra to get energy and flux arrays
    # And this will write out the csv files, so we have to remove them again...
    energy_single, flux_single = flux.write_spec_single(
        spec_single,
        Sizewell=site == "sizewell",
        Hartlepool=site == "hartlepool",
    )
    energy_multiple, flux_multiple = flux.write_spec_multiple(
        spec_multiple,
        Sizewell=site == "sizewell",
        Hartlepool=site == "hartlepool",
    )
    os.remove(f"{site.capitalize()}_single_0.5.csv")
    os.remove(f"{site.capitalize()}_multiple.csv")

    # Plot both flux spectra on one graph
    # This will save the plot to a PNG file, and doesn't return anything.
    flux.multiple_single_plot(
        energy_single,
        flux_single,
        energy_multiple,
        flux_multiple,
    )

    # Delete the output file
    # (note this will be called Sizewell regardless of the site)
    os.remove("Multiple_Single_comp_Sizewell_0.5.png")

    print("--- Multiple flux plot test passed ---")


def _test_multiple_fluxes(site="sizewell"):
    """Test plotting multiple flux spectra for different removal times."""
    print(f"--- Testing multiple fluxes for {site.capitalize()} ---")

    # Create the multiple flux spectra plots
    # This function doesn't actually plot, it just returns a list of spectra.
    # All the removal times are hardcoded in the function.
    sums_s = plotting.multiple_fluxes(
        Sizewell=site == "sizewell",
        Hartlepool=site == "hartlepool",
    )

    # Check the returned list of spectra
    assert isinstance(sums_s, ROOT.TList), "Returned object is not a TList"
    assert len(sums_s) == 5, "Returned list has wrong length"
    assert all(isinstance(spec, ROOT.TH1D) for spec in sums_s), (
        "Not all items in list are TH1D"
    )

    # Now plot the multiple spectra
    # This will pop up the ROOT window and save the plot to a pdf, but doesn't return anything.
    plotting.plot_multiple(sums_s)

    # Delete the output file
    # (again this will be called Sizewell regardless of the site)
    os.remove("Sizewell_MultipleCasks.pdf")

    # Calculate the flux at 40m for each spectrum
    if site == "sizewell":
        fluxes = [
            13067897209.144945,  # NOTE this one is the same as in _test_flux_calc
            6381971926.645747,
            1168957368.2360244,
            836003017.479554,
            653307209.876354,
        ]
    elif site == "hartlepool":
        fluxes = [
            1443268305.2325196,  # NOTE this one is the same as in _test_flux_calc
            1072427745.6884356,
            723865499.2550066,
            630534707.9926968,
            495842436.63382983,
        ]
    for spec, flux_ref in zip(sums_s, fluxes):
        flux_40 = flux.flux_calc(spec, distance_m=40)
        assert isinstance(flux_40, float), "Flux is not a float"
        assert flux_40 == flux_ref, "Flux value does not match"

    print("--- Multiple fluxes test passed ---")


def _test_sampling(site="sizewell", removal_times=[0.5, 1, 5, 10, 20]):
    """Test plotting a sample spectrum for given removal times."""
    print(f"--- Testing sampling for {site.capitalize()} ---")

    # Create the multiple cask spectrum
    if site == "sizewell":
        spec_multiple = plotting.plot_multiple_casks_sizewell(removal_times=None)
    elif site == "hartlepool":
        spec_multiple = plotting.plot_multiple_casks_hartlepool(removal_times=None)

    # Plot a sample of the spectrum
    # This will pop up the ROOT window and save the plot to a pdf, but doesn't return anything.
    # It will also save the spectrum to a CSV file, but again it's not really a proper CSV
    # just a list of randomly sampled values.
    plotting.plot_sample(spec_multiple)

    # Open the sample file and check the contents
    # Note 1 million samples is hardcoded in plot_sample
    with open("sampled_spectrum.csv") as f:
        lines = f.readlines()
        samples = [float(line.strip()) for line in lines]
    assert len(samples) == 1000000, "Wrong number of samples in CSV"

    # Compare to reference file
    # sample.sample() uses GetRandom(), but the output seems to be deterministic.
    # So we can compare to a reference file.
    with open(
        f"../../tests/test_data/{site.capitalize()}_sampled_spectrum.csv",
    ) as f:
        lines = f.readlines()
        samples_ref = [float(line.strip()) for line in lines]
    assert len(samples_ref) == 1000000, "Wrong number of samples in reference CSV"
    assert samples == samples_ref, "Sampled spectrum does not match reference"

    # Delete the output file (again the plot will be called Sizewell regardless of the site)
    os.remove("sampled_spectrum.csv")
    os.remove("Sizewell_Sampled.pdf")

    print("--- Sampling test passed ---")


def test_sizewell_commandline():
    """Run all tests for Sizewell site."""
    _test_single_cask(site="sizewell", removal_times=[0.5, 1, 5, 10, 20])
    _test_multiple_casks(site="sizewell")
    _test_flux_calc(site="sizewell", removal_times=[0.5, 1, 5, 10, 20])
    _test_multiple_plot(site="sizewell", removal_times=[0.5, 1, 5, 10, 20])
    _test_multiple_fluxes(site="sizewell")
    _test_sampling(site="sizewell", removal_times=[0.5, 1, 5, 10, 20])


def test_hartlepool_commandline():
    """Run all tests for Hartlepool site."""
    _test_single_cask(site="hartlepool", removal_times=[0.5, 1, 5, 10, 20])
    _test_multiple_casks(site="hartlepool")
    _test_flux_calc(site="hartlepool", removal_times=[0.5, 1, 5, 10, 20])
    _test_multiple_plot(site="hartlepool", removal_times=[0.5, 1, 5, 10, 20])
    _test_multiple_fluxes(site="hartlepool")
    _test_sampling(site="hartlepool", removal_times=[0.5, 1, 5, 10, 20])


if __name__ == "__main__":
    # Run the tests using the same inputs as the command line script
    test_sizewell_commandline()

    # Now run for Hartlepool as well
    test_hartlepool_commandline()

    print()
    print("All tests passed!")
