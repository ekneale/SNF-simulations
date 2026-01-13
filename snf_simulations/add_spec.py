import ROOT

import numpy as np


def add_spec(spectra):
    # adding spectra and setting limit to beyond the highest antineutrino energy in the database

    # TODO: removed for now: it doesn't seem to work, and it would break the tests.
    # for i in range(len(spectra)):
    #     spectra[i].GetXaxis().SetLimits(0, 6000)

    hsum = spectra[0].Clone("combined")
    hsum.Reset()

    hsum.Merge(spectra)
    hsum.SetTitle("Total Spectrum")

    return hsum
