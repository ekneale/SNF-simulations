import numpy as np
import ROOT


from snf_simulations import scale, add_spec, load_spec, load_and_scale, load_antineutrino_data

def AddDecays(t, parent_prop0, parent_half_life, daughter_half_life, total_mass = 1000, BR=1):

    #calculation of additional contributing isotopes that have been created from the decay of their parent isotope

    m_parent_initial = total_mass * parent_prop0 

    lambdaA = np.log(2) / parent_half_life
    lambdaB = np.log(2) / daughter_half_life

    m_daughter_new = BR*(lambdaA / (lambdaB - lambdaA)) * m_parent_initial *(np.exp(-lambdaA*t) - np.exp(-lambdaB*t))

    return m_daughter_new



def TotSpec(cask_name,removal_time =0,total_m =1000, max_E=6000, Sr90_prop=0, Y90_prop=0, Pu241_prop=0,Cs137_prop=0, Am242_prop=0, Cs135_prop=0, I129_prop=0, Np239_prop=0,Tc99_prop=0, Zr93_prop=0, Ce144_prop=0, Kr88_prop=0, Pr144_prop=0, Rb88_prop=0, Rh106_prop=0, Ru106_prop=0):

    #total mass is in kg
    #mass is in kg
    #all times given in years

    Sr90, Y90,Pu241, Cs137, Am242, Cs135, I129, Np239, Tc99, Zr93, Ce144, Kr88, Pr144, Rb88, Rh106, Ru106 = load_antineutrino_data.load_antineutrino_data()

    #name of isotopes in database

    isotopes = ["Sr90", "Y90","Pu241", "Cs137", "Am242", "Cs135", "I129", "Np239", "Tc99", "Zr93", "Ce144", "Kr88", "Pr144", "Rb88", "Rh106", "Ru106"]

    #inputted proportions of isotopes defined by the user

    proportions = [Sr90_prop, Y90_prop, Pu241_prop,Cs137_prop, Am242_prop, Cs135_prop, I129_prop,Np239_prop, Tc99_prop, Zr93_prop, Ce144_prop, Kr88_prop, Pr144_prop, Rb88_prop, Rh106_prop, Ru106_prop]

    #relative atomic masses

    mr = [90,90,241,137,242,135,129,239,99,93, 144,88,144,88,106,106]

    #half lives of isotopes

    half_life_yrs = [28.91, 0.0073, 14.329, 30.08, 0.0018, 2.4e6,1.57e7, 0.00645, 2.11e5,1.53e6, 0.781, 0.000322, 3.29e-5, 3.381e-5, 9.513e-7, 1.01]

    #calculated mass of isotopes from inputted proportions

    mass = []
    for i in range(len(proportions)):
        mass.append(proportions[i] * total_m)

    spectra = ROOT.TList()
    
    #adding the scaled spectra of each isotope to the list that will later be added together

    spectra.Add(load_and_scale.load_equal_scaled(Sr90, 546,(isotopes[0] + str(removal_time)+ str(cask_name)), isotopes[0], mass[0], mr[0], half_life_yrs[0], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Y90, 2278,(isotopes[1] + str(removal_time)+ str(cask_name)), isotopes[1], mass[1], mr[1], half_life_yrs[1], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Pu241, 20,(isotopes[2] + str(removal_time)+ str(cask_name)), isotopes[2], mass[2], mr[2], half_life_yrs[2], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Cs137, 1175,(isotopes[3] + str(removal_time)+ str(cask_name)), isotopes[3], mass[3], mr[3], half_life_yrs[3], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Am242, 664,(isotopes[4] + str(removal_time)+ str(cask_name)), isotopes[4], mass[4], mr[4], half_life_yrs[4], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Cs135, 268,(isotopes[5] + str(removal_time)+ str(cask_name)), isotopes[5], mass[5], mr[5], half_life_yrs[5], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(I129, 149,(isotopes[6] + str(removal_time)+ str(cask_name)), isotopes[6], mass[6], mr[6], half_life_yrs[6], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Np239, 714,(isotopes[7] + str(removal_time)+ str(cask_name)), isotopes[7], mass[7], mr[7], half_life_yrs[7], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Tc99, 440,(isotopes[8] + str(removal_time)+ str(cask_name)), isotopes[8], mass[8], mr[8], half_life_yrs[8], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Zr93, 1000,(isotopes[9] + str(removal_time)+ str(cask_name)), isotopes[9], mass[9], mr[9], half_life_yrs[9], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Ce144, 318,(isotopes[10] + str(removal_time)+ str(cask_name)), isotopes[10], mass[10], mr[10], half_life_yrs[10], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Kr88, 2918,(isotopes[11] + str(removal_time)+ str(cask_name)), isotopes[11], mass[11], mr[11], half_life_yrs[11], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Pr144, 2997,(isotopes[12] + str(removal_time)+ str(cask_name)), isotopes[12], mass[12], mr[12], half_life_yrs[12], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Rb88, 5301,(isotopes[13] + str(removal_time)+ str(cask_name)), isotopes[13], mass[13], mr[13], half_life_yrs[13], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Rh106, 3541,(isotopes[14] + str(removal_time)+ str(cask_name)), isotopes[14], mass[14], mr[14], half_life_yrs[14], removal_time))
    spectra.Add(load_and_scale.load_equal_scaled(Ru106, 39,(isotopes[15] + str(removal_time)+ str(cask_name)), isotopes[15], mass[15], mr[15], half_life_yrs[15], removal_time))

    #adding on the spectrum of any newly created isotopes from previous decays
    #all of these decay chains have a branching ratio of 1
    #if any additional isotopes were to be added with decay chains involving more beta emitting isotopes then they can be added here

    if removal_time != 0:
        extra_Y90 = AddDecays(t = removal_time, parent_prop0=Sr90_prop, parent_half_life=half_life_yrs[0], daughter_half_life=half_life_yrs[1], total_mass=total_m)
        spectra.Add(load_and_scale.load_equal_scaled(Y90, 2278,"additional Y90" + str(removal_time)+ str(cask_name), isotopes[1], extra_Y90, mr[1], half_life_yrs[1], 0 ))

        extra_Pr144 = AddDecays(t = removal_time, parent_prop0=Ce144_prop, parent_half_life=half_life_yrs[10], daughter_half_life=half_life_yrs[12], total_mass=total_m)
        spectra.Add(load_and_scale.load_equal_scaled(Pr144, 2997, "additional Pr144" + str(removal_time)+ str(cask_name), isotopes[12], extra_Pr144, mr[12], half_life_yrs[12], 0 ))

        extra_Rb88 = AddDecays(t = removal_time, parent_prop0=Kr88_prop, parent_half_life=half_life_yrs[11], daughter_half_life=half_life_yrs[13], total_mass=total_m)
        spectra.Add(load_and_scale.load_equal_scaled(Rb88, 5301, "additional Rb88" + str(removal_time)+ str(cask_name), isotopes[13], extra_Rb88, mr[13], half_life_yrs[13], 0 ))    

        extra_Rh106 = AddDecays(t = removal_time, parent_prop0=Ru106_prop, parent_half_life=half_life_yrs[15], daughter_half_life=half_life_yrs[14], total_mass=total_m)
        spectra.Add(load_and_scale.load_equal_scaled(Rh106, 3541, "additional Rh106" + str(removal_time)+ str(cask_name), isotopes[14], extra_Rh106, mr[14], half_life_yrs[14], 0 ))
   
    total_spec = add_spec.add_spec(spectra)
    total_spec.GetXaxis().SetRangeUser(0,max_E)
    return total_spec




    
