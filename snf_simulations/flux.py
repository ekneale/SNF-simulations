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

def write_spec(total_spec, Hartlepool = False, Sizewell = False):
    if Hartlepool == True:
        reactor = "Hartlepool"
    if Sizewell == True:
        reactor = "Sizewell"

    energy = [] # open array for energy
    flux = []
    n_bins = total_spec.GetNbinsX() #total no of bins 
    #iterating through bins assigning the centre as the energy and content as the flux 
    for i in range(1, n_bins +1):
        energy1 = total_spec.GetBinCenter(i)
        flux1 = total_spec.GetBinContent(i)
        energy.append(energy1)
        flux.append(flux1)
        #saving energy and flux data as a csv file 
# change name fir what reactor and cooling time have been chosen 
    with open(f"{reactor}_multiple_20.csv", mode= 'w', newline='') as file:
        
        file.write("\"energy\": "+str(energy)+",\n")
        file.write("\"(flux)\": "+str(flux)+",\n")
        print("Energy and Flux saved to csv")
    