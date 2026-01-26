from array import array

import numpy as np
import ROOT


def calculate_flux(spec, distance):
    """Calculate the antineutrino flux at a given distance from the source.

    Args:
        spec (ROOT.TH1D): The total antineutrino spectrum histogram.
        distance (float): Distance from the source in meters.

    Returns:
        flux (float): Antineutrino flux at the given distance in cm^{-2} s^{-1}.

    """
    # The IBD reaction threshold is 1800 keV, so we integrate above this energy.
    # This gives the total flux of antineutrinos per second above the threshold.
    total_flux = spec.Integral(1801, 6000)

    # Calculate the flux at the given distance, assuming it's a point source emitting
    # isotropically in all directions.
    # Note we convert to cm, to get the flux is in cm^-2 s^-1
    flux = (1 / (4 * np.pi * (distance * 100) ** 2)) * (total_flux)
    return flux


def calculate_event_rate(
    flux,
    lower_efficiency=0.2,
    upper_efficiency=0.4,
):
    """Calculate the expected event rate in the VIDARR detector for a given flux.

    Args:
        flux (float): Antineutrino flux in cm^-2 s^-1.
        lower_efficiency (float): Lower limit on detection efficiency. Default is 0.2.
        upper_efficiency (float): Upper limit on detection efficiency. Default is 0.4

    Returns:
        rate_lower, rate_upper (tuple of floats): event rates in s^{-1}
            for the lower and upper efficiency limits.

    """
    # Calculate the number of target protons in the detector
    # VIDARR is a plastic scintillator detector with a volume of (1.52 x 1.52 x 0.7) m^3
    detector_volume = 1.52 * 1.52 * 0.7  # m^3
    detector_volume = detector_volume * 1e6  # convert to cm^3
    proton_density = 4.6e22  # number density of protons in cm^-3
    number_of_protons = detector_volume * proton_density

    # Calculate the event rate using the flux and number of target protons
    cross_section = 1e-44  # cm^2, approximate IBD cross-section
    event_rate = number_of_protons * cross_section * flux

    # Apply detection efficiency
    rate_lower = event_rate * lower_efficiency
    rate_upper = event_rate * upper_efficiency
    return rate_lower, rate_upper


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
