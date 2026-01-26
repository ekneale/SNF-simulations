from array import array

import numpy as np
import ROOT


def flux_calc(
    total_spec,
    distance_m,
):
    total_flux_above_threshold = total_spec.Integral(
        1801,
        6000,
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
