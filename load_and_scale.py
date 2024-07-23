import numpy as np
import ROOT

import scale
import load_spec

def load_equal_scaled(data,max_E,named, isotope, m,mr,half_life_yrs,removal_time, min_E=0):
    spec = load_spec.load_equal(named,isotope, data[:,0], data[:,1], max_E, min_E)
    spec_scaled = scale.scale(spec, m, mr, half_life_yrs, removal_time)
    return spec_scaled

