from define_proportions import TotSpec, AddDecays
import numpy as np
from sample import sample
import ROOT
import load_spec
import load_and_scale

#total1 = TotSpec(removal_time=0, I129_prop=259e-6, Cs135_prop=1300e-6, Cs137_prop=1300e-6, Sr90_prop=794e-6, Zr93_prop=3639e-6, Tc99_prop=799e-6, Am242_prop=484e-6,Ce144_prop=2469e-6, Pr144_prop=1161e-6, Rh106_prop=484e-6, Ru106_prop=2404e-6)

total10 = TotSpec(removal_time=10, I129_prop=259e-6, Cs135_prop=1300e-6, Cs137_prop=1300e-6, Sr90_prop=794e-6, Zr93_prop=3639e-6, Tc99_prop=799e-6, Am242_prop=484e-6,Ce144_prop=2469e-6, Pr144_prop=1161e-6, Rh106_prop=484e-6, Ru106_prop=2404e-6)

total30 = TotSpec(removal_time=30, I129_prop=259e-6, Cs135_prop=1300e-6, Cs137_prop=1300e-6, Sr90_prop=794e-6, Zr93_prop=3639e-6, Tc99_prop=799e-6, Am242_prop=484e-6,Ce144_prop=2469e-6, Pr144_prop=1161e-6, Rh106_prop=484e-6, Ru106_prop=2404e-6)

total30.SetLineColor(ROOT.kGreen)

total10.SetLineColor(ROOT.kRed)

#total1.Draw("hist")

total10.Draw("hist")

total30.Draw("hist same")




input("exit")