import ROOT

import numpy as np

def extend_spec(NewMax, NewBins, og_spec):
    extended = ROOT.TH1D("extended" +str(NewMax), "", NewBins, 0, NewMax)

    for i in range(og_spec.GetNbinsX()):
        bin_center = og_spec.GetBinCenter(i)
        content = og_spec.GetBinContent(i)
        error = og_spec.GetBinError(i)
        new_bin = og_spec.FindBin(bin_center)
        extended.SetBinContent(new_bin, content)
        extended.SetBinError(new_bin, error)

    return extended

def add_byStack(spectra):
    maxima = []
    for i in range(len(spectra)):
        maxima.append(spectra[i].GetXaxis().GetXmax())
        

    for i in range(len(spectra)):
        spectra[i].GetXaxis().SetLimits(0, max(maxima))

    hsum = ROOT.THStack("spectra","")

    for i in range(len(spectra)):
        hsum.Add(spectra[i])

    return hsum

def add_byMerge(spectra):

    maxima = []
    for i in range(len(spectra)):
        maxima.append(spectra[i].GetXaxis().GetXmax())
        

    for i in range(len(spectra)):
        spectra[i].GetXaxis().SetLimits(0, max(maxima))

    hsum = spectra[14].Clone("combined")
    hsum.Reset()

    hsum.Merge(spectra)
    hsum.SetTitle("Total Spectrum")
    hsum.GetXaxis().SetTitle("Energy [keV]")
    hsum.GetYaxis().SetTitle("Relative Flux [s$^{-1}$]")
 
    return hsum


    
    
