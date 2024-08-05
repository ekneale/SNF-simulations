import numpy as np
from snf_simulations import sample
import ROOT



def flux_calc(total_spec,distance_m,):

    total_flux_above_threshold = total_spec.Integral(1801, 6000) #total flux per second above threshold

    print(total_flux_above_threshold, "per second")

    flux = (1/(4*np.pi*(distance_m*100)**2)) * (total_flux_above_threshold)

    print(flux, "cm^-2 s^-1")
    print(flux * 60*60*24, "cm^-2 days^-1")

    N_p = (1.52 * 1.52 * 0.7) * 1e6 * 4.6*1e23

    #calulcating event rate for lower limit and upper limit on efficiency
    Rate_lower = N_p * 1e-44 * flux * 0.2

    Rate_upper = N_p * 1e-44 * flux * 0.4

    print(Rate_lower, Rate_upper, "per s")
    print(Rate_lower * 60* 60 * 24, Rate_upper * 60**2 * 24,"per day")


    return flux
