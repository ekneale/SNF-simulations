"""Command line script to run SNF simulations and generate plots."""

from array import array

import numpy as np
import ROOT

from snf_simulations.cask import get_total_spec
from snf_simulations.data import load_reactor_data
from snf_simulations.physics import calculate_event_rate, calculate_flux
from snf_simulations.sample import sample_spec
from snf_simulations.spec import add_spec, write_spec

ROOT.TH1.AddDirectory(False)  # Prevent ROOT from keeping histograms in memory

# Suppress warnings from ruff
# ruff: noqa: T201  # prints


def _get_spectra(cask_mass=10000, removal_times=None, reactor="sizewell"):
    """Create spectra for a single dry cask of fuel at multiple removal times.

    Args:
        cask_mass (float): Total mass of SNF in the cask (kg).
            Default is 10 tonnes (10,000 kg).
        removal_times (list or None): List of removal times to plot (years).
            If None, only plot time 0.
        reactor (str): Reactor name. Default is "sizewell".

    Returns:
        ROOT.TH1D: Spectrum for the last removal time in the list.

    """
    if removal_times is None:
        removal_times = [0]

    # Load the isotope proportions for the specified reactor
    # This will raise a ValueError if the reactor name is invalid
    proportions = load_reactor_data(reactor)

    # Create specs
    spectra = ROOT.TList()
    for removal_time in removal_times:
        spec = get_total_spec(
            cask_name=f"{reactor}_cask_{cask_mass}_{removal_time}",
            isotope_proportions=proportions,
            total_mass=cask_mass,
            removal_time=removal_time,
        )
        spectra.Add(spec)

    return spectra


def main(reactor="sizewell", cask_mass=10000):
    """Run SNF simulations and generate plots."""
    ################################################
    # 1) Plot the spectrum for a single dry cask of fuel at different removal times.
    #    One 10 tonne cask, removal times of 0, 0.5, 1, 5, 10, 20 years.
    print("Generating single cask spectra at different cooling times...")
    removal_times = [0, 0.5, 1, 5, 10, 20]

    # Create the spectra
    spectra = _get_spectra(
        cask_mass=cask_mass,
        removal_times=removal_times,
        reactor=reactor,
    )

    # Calculate and print flux and event rates for the 0.5 year removal time
    spec_single_05 = spectra[1]
    print()
    print(f"Single cask flux at 40 m for {reactor.capitalize()} after 0.5 years:")
    flux = calculate_flux(spec_single_05, 40)
    print(f"Flux: {flux:.3e} cm^-2 s^-1")
    print(f"Flux: {flux * 60 * 60 * 24:.3e} cm^-2 days^-1")
    rate_lower, rate_upper = calculate_event_rate(flux, 0.2, 0.4)
    print(f"Event rate: {rate_lower:.3e} to {rate_upper:.3e} per s")
    print(
        f"Event rate: {rate_lower * 60 * 60 * 24:.3f} to {rate_upper * 60**2 * 24:.3f} per day",
    )

    # Save 0.5 year spectrum to CSV
    print()
    print("Writing single cask spectrum data to CSV...")
    filename = f"{reactor}_single.csv"
    data_single_05 = write_spec(spec_single_05, output_filename=filename)
    print(f"Saved to {filename}")

    # Create ROOT canvas and legend
    c = ROOT.TCanvas("single_cask", "Total Spectrum", 1200, 600)
    legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)

    # Add spectra to canvas
    colors = [
        ROOT.kBlue,
        ROOT.kRed,
        ROOT.kGreen,
        ROOT.kBlack,
        ROOT.kViolet,
        ROOT.kOrange,
    ]
    for i in range(len(spectra)):
        spectra[i].SetLineColor(colors[i % len(colors)])
        legend.AddEntry(
            spectra[i],
            f"{removal_times[i]} years since removal from core",
        )
        spectra[i].Draw("hist same")
    legend.Draw()

    # Format plot
    c.SetLogy()
    spectra[0].SetTitle("")
    spectra[0].GetXaxis().SetTitle("Energy [keV]")
    spectra[0].GetXaxis().SetLabelSize(0.05)
    spectra[0].GetXaxis().SetTitleSize(0.05)
    spectra[0].GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    spectra[0].GetYaxis().SetLabelSize(0.05)
    spectra[0].GetYaxis().SetTitleSize(0.05)

    # Display the plot
    c.Update()
    input("Displaying plot...")  # pause to view plot

    # Save the plot as a PDF
    c.SaveAs(f"{reactor.capitalize()}_single.pdf")

    ################################################
    # 2) Calculate total spectra for multiple casks
    #    40 casks in total, each 10 tonnes, split between 4 different removal times.
    #    Since we're combining the spectra at the end, instead of ten 10-tonne casks
    #    at each removal time, we just simulate one 100-tonne cask.
    print()
    print("Generating multiple cask spectra...")
    casks_per_removal = 10
    if reactor == "sizewell":
        removal_times = [0.5, 5, 10, 20]
    elif reactor == "hartlepool":
        removal_times = [3, 7, 15, 19]

    # Create the spectra and combine them
    spectra = _get_spectra(
        cask_mass=cask_mass * casks_per_removal,
        removal_times=removal_times,
        reactor=reactor,
    )
    spec_multiple = add_spec(spectra)
    spec_multiple.SetTitle("Total spectrum for all casks")

    # Save to CSV
    print()
    print("Writing single cask spectrum data to CSV...")
    filename = f"{reactor}_multiple.csv"
    data_multiple = write_spec(spec_multiple, output_filename=filename)
    print(f"Saved to {filename}")

    # Calculate and print flux and event rates for the combined spectrum
    print()
    print(f"Multiple cask flux at 40 m for {reactor.capitalize()}:")
    flux = calculate_flux(spec_multiple, 40)
    print(f"Flux: {flux:.3e} cm^-2 s^-1")
    print(f"Flux: {flux * 60 * 60 * 24:.3e} cm^-2 days^-1")
    rate_lower, rate_upper = calculate_event_rate(flux, 0.2, 0.4)
    print(f"Event rate: {rate_lower:.3e} to {rate_upper:.3e} per s")
    print(
        f"Event rate: {rate_lower * 60 * 60 * 24:.3f}"
        f" to {rate_upper * 60**2 * 24:.3f} per day",
    )

    ################################################
    # 3) Plot both flux spectra on one graph
    # Create ROOT canvas and legend
    c = ROOT.TCanvas("multiple_single", "true_neutrino_energy", 1200, 600)
    legend = ROOT.TLegend(0.7, 0.15, 0.9, 0.3)

    # Add data to the canvas
    energy_single, flux_single = data_single_05[:, 0], data_single_05[:, 1]
    graph_single = ROOT.TGraph(
        len(energy_single),
        array("d", energy_single),
        array("d", flux_single),
    )
    graph_single.SetLineColor(ROOT.kBlue)
    graph_single.Draw("APL")
    legend.AddEntry(graph_single, f"{reactor.capitalize()}, single cask (0.5 yrs)")

    energy_multiple, flux_multiple = data_multiple[:, 0], data_multiple[:, 1]
    graph_multiple = ROOT.TGraph(
        len(energy_multiple),
        array("d", energy_multiple),
        array("d", flux_multiple),
    )
    graph_multiple.SetLineColor(ROOT.kRed)
    graph_multiple.Draw("same L")
    legend.AddEntry(graph_multiple, f"{reactor.capitalize()}, multiple casks")
    legend.Draw()

    # Format plot
    c.SetLogy()
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # Display the plot
    c.Update()
    input("Displaying plot...")  # pause to view plot

    # Save the plot as a PDF
    c.SaveAs(f"Multiple_Single_comp_{reactor.capitalize()}_0.5.png")

    ################################################
    # 4) Plot multiple flux spectra for different cooling times.
    # As in part 2, but now we add extra cooling times to see how the spectrum evolves.
    print()
    print("Generating multiple cask spectra for different cooling times...")
    cask_mass = 10000
    casks_per_removal = 10
    if reactor == "sizewell":
        removal_times = [0.5, 5, 10, 20]
    elif reactor == "hartlepool":
        removal_times = [3, 7, 15, 19]
    cooling_times = [0, 1, 5, 10, 20]

    # Create the combined spectra from all casks for each cooling time
    spectra = ROOT.TList()
    for cooling_time in cooling_times:
        time_spectra = _get_spectra(
            cask_mass=cask_mass * casks_per_removal,
            removal_times=np.array(removal_times) + cooling_time,
            reactor=reactor,
        )
        total_spec = add_spec(time_spectra)
        spectra.Add(total_spec)

    # Calculate and print flux and event rates for each cooling time
    print()
    for i in range(len(spectra)):
        print(
            f"Multiple cask flux at 40 m for {reactor.capitalize()}"
            f" after {cooling_times[i]} years:",
        )
        flux = calculate_flux(spectra[i], 40)
        print(f"Flux: {flux:.3e} cm^-2 s^-1")
        print(f"Flux: {flux * 60 * 60 * 24:.3e} cm^-2 days^-1")
        rate_lower, rate_upper = calculate_event_rate(flux, 0.2, 0.4)
        print(f"Event rate: {rate_lower:.3e} to {rate_upper:.3e} per s")
        print(
            f"Event rate: {rate_lower * 60 * 60 * 24:.3f}"
            f" to {rate_upper * 60**2 * 24:.3f} per day",
        )
        print()

    # Create ROOT canvas and legend
    c = ROOT.TCanvas("multiple", "Total Spectrum", 1200, 600)
    legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)

    # Add data to the canvas
    colors = [ROOT.kBlue, ROOT.kGreen, ROOT.kViolet, ROOT.kBlack, ROOT.kRed]
    for i in range(len(spectra)):
        spectra[i].SetLineColor(colors[i])
        spectra[i].Draw("hist same")
        legend.AddEntry(
            spectra[i],
            "Spectrum after " + str(cooling_times[i]) + " years",
        )
    legend.Draw()

    # Format plot
    c.SetLogy()
    spectra[0].SetTitle("")
    spectra[0].GetXaxis().SetTitle("Energy [keV]")
    spectra[0].GetXaxis().SetLabelSize(0.05)
    spectra[0].GetXaxis().SetTitleSize(0.05)
    spectra[0].GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    spectra[0].GetYaxis().SetLabelSize(0.05)
    spectra[0].GetXaxis().SetTitleSize(0.05)

    # Display the plot
    c.Update()
    input("Displaying plot...")  # pause to view plot

    # Save the plot as a PDF
    c.SaveAs(f"{reactor.capitalize()}_MultipleCasks.pdf")

    ################################################
    # 5) Sample the multiple cask spectrum to simulate detector observations
    print()
    print("Sampling multiple cask spectrum...")

    # Take 1 million samples from the multiple cask spectrum, and save to CSV
    sampled = sample_spec(
        spec_multiple,
        samples=1000000,
        output_filename=f"{reactor.capitalize()}_sampled_spectrum.csv",
    )

    # Create ROOT canvas and legend
    c = ROOT.TCanvas("samples", "Total Spectrum", 1200, 600)
    c.SetLogy()

    # Plot a histogram of the sampled data
    h = ROOT.TH1D("sample", "", 6000, 0, 6000)
    for i in range(1000000):
        h.Fill(sampled[i])
    h.Draw("hist")

    # Format plot
    h.SetStats(0)
    h.SetTitle("")
    h.GetXaxis().SetTitle("Energy [keV]")
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    h.GetYaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)

    # Display the plot
    c.Update()
    input("Displaying plot...")  # pause to view plot

    # Save the plot as a PDF
    c.SaveAs(f"{reactor.capitalize()}_Sampled.pdf")


if __name__ == "__main__":
    main()
