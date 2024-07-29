import numpy as np

import ROOT

def add_casks(casks):
    cask_sum = casks[0].Clone("casks")
    cask_sum.Reset()

    cask_sum.Merge(casks)
    cask_sum.SetTitle("Total spectrum for all casks")

    c= ROOT.TCanvas

    cask_sum.SetTitle("Resultant Spectrum of Several Dry Casks")
    cask_sum.GetXaxis.SetTitle("Energy [keV]")
    cask_sum.GetYaxis.SetTitle("Relative Flux [s #^{-1}]")
    
    return cask_sum
