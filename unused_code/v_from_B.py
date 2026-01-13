import ROOT
import numpy as np


def v_from_B(Energy, dN, isotope):
    c = ROOT.TCanvas()

    n_bins = len(Energy) - 1

    bin_width = []

    for i in range(n_bins):
        bin_width.append(Energy[i + 1] - Energy[i])

    new_widths = np.array(bin_width)[::-1]

    edges = [0.0]

    for i in range(n_bins):
        edges.append(edges[i] + new_widths[i])

    hist = ROOT.TH1D(isotope, isotope, n_bins, np.array(edges))

    for i in range(len(Energy)):
        hist.Fill(max(Energy) - Energy[i], dN[i])

    hist.SetStats(0)

    hist.GetXaxis().SetTitle("Energy [keV]")
    hist.SetTItle("Antineutrino Energy Spectrum")
    c.Update()

    return hist
