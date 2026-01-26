from array import array

import ROOT

from .cask import get_total_spec
from .data import load_reactor_data
from .sample import sample
from .spec import add_spec

# the things that have been changed so far:
# for the single casks, the energy and flux will save to a csv file. You can pick which cooling time by changing the number in the spectra.At().
# You can pick the reactor in command line, the name of the reactor is now a variable so you don't have to manually change the pdf name .
# The multiple flux function has been updated so that sizewell is defined as well as Hartlepool. The plot spectrum function also has the reactor variable.
# Two seperate functions are defined for plotting multiple casks for either sizewell or hartlepool as they have different removal times so energy/flux csv can be found per time
# At the moment you have to manually change the csv file name for what reactor, cooling time etc but eventually this would be done automatically.


def plot_single_cask(reactor, removal_times):
    """Plot single cask spectra for given reactor and removal times.

    Args:
        reactor (str): Reactor name, either "sizewell" or "hartlepool".
        removal_times (list): List of removal times in years.

    Returns:
        ROOT.TH1D: The first spectrum in the list for further use.

    """
    # Load the isotope proportions for the specified reactor
    # This will raise a ValueError if the reactor name is invalid
    proportions = load_reactor_data(reactor)

    # Create canvas and legend
    c = ROOT.TCanvas("c", "Total Spectrum", 1200, 600)
    c.SetLogy()
    legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)

    # Get initial spectrum at removal time 0 for reference
    initial = get_total_spec(
        cask_name="main",
        isotope_proportions=proportions,
        total_mass=100000,
        removal_time=0,
    )
    initial.SetLineColor(ROOT.kBlue)
    legend.AddEntry(initial, "1 Day since removal from core")
    initial.Draw("hist")

    # Get spectra for specified removal times
    spectra = ROOT.TList()
    colors = [ROOT.kRed, ROOT.kGreen, ROOT.kBlack, ROOT.kViolet, ROOT.kOrange]
    for removal_time, color in zip(removal_times, colors):
        spec = get_total_spec(
            cask_name="main",
            isotope_proportions=proportions,
            total_mass=100000,
            removal_time=removal_time,
        )
        spec.SetLineColor(color)
        spectra.Add(spec)
    for i in range(len(spectra)):
        legend.AddEntry(
            spectra[i],
            str(removal_times[i]) + " Years since removal from core",
        )
        spectra[i].Draw("hist same")

    # Format plot
    initial.SetTitle("")
    initial.GetXaxis().SetTitle("Energy [keV]")
    initial.GetXaxis().SetLabelSize(0.05)
    initial.GetXaxis().SetTitleSize(0.05)
    initial.GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    initial.GetYaxis().SetLabelSize(0.05)
    initial.GetXaxis().SetTitleSize(0.05)

    # Display and save plot
    c.Update()
    legend.Draw()
    input("exit")  # pause to view plot
    c.SaveAs(f"{reactor.capitalize()}_Spectra_0.5.pdf")  # change no for cooling time

    return spectra.At(0)  # needs .At() as its a Tlist not a python array


def plot_multiple_casks_sizewell(removal_times):
    # plotting 10 dry casks x 10tonnes  x 4 cooling times for Sizewell PWR
    proportions = load_reactor_data("sizewell")
    removal_times = [0.5, 5, 10, 20]
    casks = ROOT.TList()

    for i in range(len(removal_times)):
        casks.Add(
            get_total_spec(
                cask_name="Sizewell",
                isotope_proportions=proportions,
                total_mass=100000,
                removal_time=removal_times[i],
            ),
        )

    total_sizewell = add_spec(casks)
    total_sizewell.SetTitle("Total spectrum for all casks")

    # c.SaveAs("Sizewell_MultipleCasks.pdf")

    return total_sizewell


def plot_multiple_casks_hartlepool(removal_times):
    proportions = load_reactor_data("hartlepool")
    removal_times = [3, 7, 15, 19]
    casks = ROOT.TList()

    for i in range(len(removal_times)):
        casks.Add(
            get_total_spec(
                cask_name="Hartlepool",
                isotope_proportions=proportions,
                total_mass=100000,
                removal_time=removal_times[i],
            ),
        )

    total_hartlepool = add_spec(casks)
    total_hartlepool.SetTitle("Total spectrum for all casks")

    # c.SaveAs("Hartlepool_MultipleCasks.pdf")

    return total_hartlepool


# just using the plot_multiple_casks results in a graphics error due to the amount spectra being added so plot has to be saved separately and then plotted


def plot(spectrum, reactor):
    c = ROOT.TCanvas("c", "Total Spectrum", 1200, 600)
    c.SetLogy()

    spectrum.Draw("histS")
    spectrum.SetTitle("")
    spectrum.GetXaxis().SetTitle("Energy [keV]")
    spectrum.GetXaxis().SetLabelSize(0.05)
    spectrum.GetXaxis().SetTitleSize(0.05)
    spectrum.GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    spectrum.GetYaxis().SetLabelSize(0.05)
    spectrum.GetXaxis().SetTitleSize(0.05)

    c.Update()
    input("exit")  # pause to view plot
    c.SaveAs(f"{reactor.capitalize()}_casks.pdf")


def multiple_fluxes(reactor):
    """Plot multiple flux spectra for different cooling times.

    Args:
        reactor (str): Reactor name, either "sizewell" or "hartlepool".

    Returns:
        ROOT.TList: List of summed spectra for different cooling times.

    """
    extra_times = [1, 5, 10, 20]
    if reactor == "hartlepool":
        removal_times = [3, 7, 15, 19]
    elif reactor == "sizewell":
        removal_times = [0.5, 5, 10, 20]
    else:
        raise ValueError("Reactor must be either 'sizewell' or 'hartlepool'")
    proportions = load_reactor_data(reactor)

    casks0 = ROOT.TList()
    casks1 = ROOT.TList()
    casks5 = ROOT.TList()
    casks10 = ROOT.TList()
    casks20 = ROOT.TList()

    for i in range(len(removal_times)):
        casks0.Add(
            get_total_spec(
                cask_name=f"{reactor.capitalize()}0",
                isotope_proportions=proportions,
                total_mass=100000,
                removal_time=removal_times[i],
            ),
        )

    sum0 = add_spec(casks0)
    sum0.SetTitle("Total spectrum for all casks")

    for i in range(len(removal_times)):
        casks1.Add(
            get_total_spec(
                cask_name=f"{reactor.capitalize()}1",
                isotope_proportions=proportions,
                total_mass=100000,
                removal_time=removal_times[i] + extra_times[0],
            ),
        )

    sum1 = add_spec(casks1)
    sum1.SetTitle("Total spectrum for all casks")

    for i in range(len(removal_times)):
        casks5.Add(
            get_total_spec(
                cask_name=f"{reactor.capitalize()}2",
                isotope_proportions=proportions,
                total_mass=100000,
                removal_time=removal_times[i] + extra_times[1],
            ),
        )

    sum5 = add_spec(casks5)
    sum5.SetTitle("Total spectrum for all casks")

    for i in range(len(removal_times)):
        casks10.Add(
            get_total_spec(
                cask_name=f"{reactor.capitalize()}3",
                isotope_proportions=proportions,
                total_mass=100000,
                removal_time=removal_times[i] + extra_times[2],
            ),
        )

    sum10 = add_spec(casks10)
    sum10.SetTitle("Total spectrum for all casks")

    for i in range(len(removal_times)):
        casks20.Add(
            get_total_spec(
                cask_name=f"{reactor.capitalize()}4",
                isotope_proportions=proportions,
                total_mass=100000,
                removal_time=removal_times[i] + extra_times[3],
            ),
        )

    sum20 = add_spec(casks20)
    sum20.SetTitle("Total spectrum for all casks")

    sums = ROOT.TList()
    sums.Add(sum0)
    sums.Add(sum1)
    sums.Add(sum5)
    sums.Add(sum10)
    sums.Add(sum20)
    return sums


# plotting separately to avoid graphics errors


def plot_multiple(multiple, reactor):
    c = ROOT.TCanvas("c", "Total Spectrum", 1200, 600)
    c.SetLogy()

    legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)

    extra_times = [0, 1, 5, 10, 20]
    multiple[0].SetLineColor(ROOT.kBlue)
    multiple[0].Draw("hist")

    legend.AddEntry(multiple[0], "Total spectrum")
    colors = [ROOT.kYellow, ROOT.kGreen, ROOT.kViolet, ROOT.kBlack, ROOT.kRed]

    for i in range(1, len(multiple)):
        multiple[i].SetLineColor(colors[i])
        multiple[i].Draw("hist same")
        legend.AddEntry(multiple[i], "Spectrum after " + str(extra_times[i]) + " years")

    multiple[0].SetTitle("")
    multiple[0].GetXaxis().SetTitle("Energy [keV]")
    multiple[0].GetXaxis().SetLabelSize(0.05)
    multiple[0].GetXaxis().SetTitleSize(0.05)
    multiple[0].GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    multiple[0].GetYaxis().SetLabelSize(0.05)
    multiple[0].GetXaxis().SetTitleSize(0.05)

    legend.Draw()
    c.Update()
    input("exit")

    c.SaveAs(f"{reactor.capitalize()}_MultipleCasks.pdf")


def plot_sample(total_spec, reactor):
    h = ROOT.TH1D("sample", "", 6000, 0, 6000)
    sampled = sample(total_spec, N=1000000)

    for i in range(1000000):
        h.Fill(sampled[i])

    c = ROOT.TCanvas("c", "Total Spectrum", 1200, 600)
    c.SetLogy()

    h.Draw("hist")
    h.SetStats(0)
    h.SetTitle("")
    h.GetXaxis().SetTitle("Energy [keV]")
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    h.GetYaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)

    c.Update()
    input("exit")
    c.SaveAs(f"{reactor.capitalize()}_Sampled.pdf")


def multiple_single_plot(
    energy_single,
    flux_single,
    energy_multiple,
    flux_multiple,
    reactor,
):
    c = ROOT.TCanvas("c", "true_neutrino_energy", 1200, 600)
    c.SetLogy()
    # c.SetLogx()
    legend = ROOT.TLegend(0.7, 0.15, 0.9, 0.3)
    graph_single = ROOT.TGraph(
        len(energy_single),
        array("d", energy_single),
        array("d", flux_single),
    )
    graph_single.SetLineColor(ROOT.kBlue)

    graph_multiple = ROOT.TGraph(
        len(energy_multiple),
        array("d", energy_multiple),
        array("d", flux_multiple),
    )
    graph_multiple.SetLineColor(ROOT.kRed)
    graph_single.Draw("APL")
    graph_multiple.Draw("same L")

    legend.AddEntry(graph_single, f"Single cask {reactor.capitalize()} 0.5 yrs")
    legend.AddEntry(graph_multiple, f"{reactor.capitalize()} multiple casks")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    c.Update()
    c.SaveAs(f"Multiple_Single_comp_{reactor.capitalize()}_0.5.png")
    input("press enter to exit")
