import numpy as np
from snf_simulations import sample
import ROOT
from array import array


def flux_calc(
    total_spec,
    distance_m,
):
    total_flux_above_threshold = total_spec.Integral(
        1801, 6000
    )  # total flux per second above threshold

    print(total_flux_above_threshold, "per second")

    flux = (1 / (4 * np.pi * (distance_m * 100) ** 2)) * (total_flux_above_threshold)

    print(flux, "cm^-2 s^-1")
    print(flux * 60 * 60 * 24, "cm^-2 days^-1")

    N_p = (1.52 * 1.52 * 0.7) * 1e6 * 4.6 * 1e23

    # calulcating event rate for lower limit and upper limit on efficiency
    Rate_lower = N_p * 1e-44 * flux * 0.2

    Rate_upper = N_p * 1e-44 * flux * 0.4

    print(Rate_lower, Rate_upper, "per s")
    print(Rate_lower * 60 * 60 * 24, Rate_upper * 60**2 * 24, "per day")

    return flux


def write_spec_multiple(total_sizewell, Hartlepool=False, Sizewell=False):
    if Hartlepool == True:
        reactor = "Hartlepool"
    if Sizewell == True:
        reactor = "Sizewell"

    energy_multiple = []
    flux_multiple = []
    n_bins = total_sizewell.GetNbinsX()

    for i in range(1, n_bins + 1):
        energy1 = total_sizewell.GetBinCenter(i) / 1e3
        flux1 = total_sizewell.GetBinContent(i)
        energy_multiple.append(energy1)
        flux_multiple.append(flux1)
        # saving energy and flux data as a csv file

    with open(f"{reactor}_multiple.csv", mode="w", newline="") as file:
        file.write('"energy": ' + str(energy_multiple) + ",\n")
        file.write('"(flux)": ' + str(flux_multiple) + ",\n")
        print("Energy and Flux saved to csv")

        return energy_multiple, flux_multiple


def write_spec_single(total_spec, Hartlepool=False, Sizewell=False):
    if Hartlepool == True:
        reactor = "Hartlepool"
    if Sizewell == True:
        reactor = "Sizewell"

    energy_single = []
    flux_single = []
    n_bins = total_spec.GetNbinsX()

    for i in range(1, n_bins + 1):
        energy1 = total_spec.GetBinCenter(i) / 1e3
        flux1 = total_spec.GetBinContent(i)
        energy_single.append(energy1)
        flux_single.append(flux1)
        # saving energy and flux data as a csv file

    # change name for what cooling time has been chosen
    with open(f"{reactor}_single_0.5.csv", mode="w", newline="") as file:
        file.write('"energy": ' + str(energy_single) + ",\n")
        file.write('"(flux)": ' + str(flux_single) + ",\n")
        print("Energy and Flux saved to csv")

        return energy_single, flux_single


def multiple_single_plot(energy_single, flux_single, energy_multiple, flux_multiple):
    c = ROOT.TCanvas("c", "true_neutrino_energy", 1200, 600)
    c.SetLogy()
    # c.SetLogx()
    legend = ROOT.TLegend(0.7, 0.15, 0.9, 0.3)
    graph_single = ROOT.TGraph(
        len(energy_single), array("d", energy_single), array("d", flux_single)
    )
    graph_single.SetLineColor(ROOT.kBlue)

    graph_multiple = ROOT.TGraph(
        len(energy_multiple), array("d", energy_multiple), array("d", flux_multiple)
    )
    graph_multiple.SetLineColor(ROOT.kRed)
    graph_single.Draw("APL")
    graph_multiple.Draw("same L")

    legend.AddEntry(graph_single, "Single cask sizewell 0.5 yrs")
    legend.AddEntry(graph_multiple, "Sizewell multiple casks")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    c.Update()
    c.SaveAs("Multiple_Single_comp_Sizewell_0.5.png")
    input("press enter to exit")
