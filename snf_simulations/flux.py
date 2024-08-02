import numpy as np
from snf_simulations import sample
import ROOT



def flux_calc(total_spec,distance_m,N=1000000):
    events = np.array(sample.sample(total_spec, N))

    total_events = total_spec.Integral() #total flux per second

    flux = (1/(4*np.pi*(distance_m*100)**2)) * (total_events / len(events))

    print(flux, "cm^-2 s^-1")

    return flux


