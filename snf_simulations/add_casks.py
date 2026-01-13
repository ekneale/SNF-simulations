import numpy as np

import ROOT


def add_casks(casks):
    cask_sum = casks[0].Clone("casks")
    cask_sum.Reset()

    cask_sum.Merge(casks)
    cask_sum.SetTitle("Total spectrum for all casks")

    return cask_sum
