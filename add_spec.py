import ROOT

import numpy as np

def add_spec(spectra):

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


    