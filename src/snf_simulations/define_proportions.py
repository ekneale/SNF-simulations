import numpy as np
import ROOT

from .load_data import load_antineutrino_data, load_isotope_data
from .spec import add_spec, load_spec


def AddDecays(
    t,
    parent_prop0,
    parent_half_life,
    daughter_half_life,
    total_mass=1000,
    BR=1,
):
    # calculation of additional contributing isotopes that have been created from the decay of their parent isotope

    m_parent_initial = total_mass * parent_prop0

    lambdaA = np.log(2) / parent_half_life
    lambdaB = np.log(2) / daughter_half_life

    m_daughter_new = (
        BR
        * (lambdaA / (lambdaB - lambdaA))
        * m_parent_initial
        * (np.exp(-lambdaA * t) - np.exp(-lambdaB * t))
    )

    return m_daughter_new


def TotSpec(
    cask_name,
    removal_time=0,
    total_m=1000,
    max_E=6000,
    Sr90_prop=0,
    Y90_prop=0,
    Pu241_prop=0,
    Cs137_prop=0,
    Am242_prop=0,
    Cs135_prop=0,
    I129_prop=0,
    Np239_prop=0,
    Tc99_prop=0,
    Zr93_prop=0,
    Ce144_prop=0,
    Kr88_prop=0,
    Pr144_prop=0,
    Rb88_prop=0,
    Rh106_prop=0,
    Ru106_prop=0,
):
    # total mass is in kg
    # mass is in kg
    # all times given in years

    # inputted proportions of isotopes defined by the user
    proportions = {
        "Sr90": Sr90_prop,
        "Y90": Y90_prop,
        "Pu241": Pu241_prop,
        "Cs137": Cs137_prop,
        "Am242": Am242_prop,
        "Cs135": Cs135_prop,
        "I129": I129_prop,
        "Np239": Np239_prop,
        "Tc99": Tc99_prop,
        "Zr93": Zr93_prop,
        "Ce144": Ce144_prop,
        "Kr88": Kr88_prop,
        "Pr144": Pr144_prop,
        "Rb88": Rb88_prop,
        "Rh106": Rh106_prop,
        "Ru106": Ru106_prop,
    }

    # Load the isotope data dicts
    isotopes = list(proportions.keys())
    atomic_masses, half_lives = load_isotope_data(isotopes)
    isotope_data = load_antineutrino_data(isotopes)

    # calculated mass of isotopes from inputted proportions
    masses = {isotope: prop * total_m for isotope, prop in proportions.items()}

    # adding the scaled spectra of each isotope to the list that will later be added together
    spectra = ROOT.TList()
    spectra.Add(
        load_spec(
            isotope_data["Sr90"],
            ("Sr90" + str(removal_time) + str(cask_name)),
            masses["Sr90"],
            atomic_masses["Sr90"],
            half_lives["Sr90"],
            removal_time,
            546,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Y90"],
            ("Y90" + str(removal_time) + str(cask_name)),
            masses["Y90"],
            atomic_masses["Y90"],
            half_lives["Y90"],
            removal_time,
            2278,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Pu241"],
            ("Pu241" + str(removal_time) + str(cask_name)),
            masses["Pu241"],
            atomic_masses["Pu241"],
            half_lives["Pu241"],
            removal_time,
            20,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Cs137"],
            ("Cs137" + str(removal_time) + str(cask_name)),
            masses["Cs137"],
            atomic_masses["Cs137"],
            half_lives["Cs137"],
            removal_time,
            1175,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Am242"],
            ("Am242" + str(removal_time) + str(cask_name)),
            masses["Am242"],
            atomic_masses["Am242"],
            half_lives["Am242"],
            removal_time,
            664,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Cs135"],
            ("Cs135" + str(removal_time) + str(cask_name)),
            masses["Cs135"],
            atomic_masses["Cs135"],
            half_lives["Cs135"],
            removal_time,
            268,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["I129"],
            ("I129" + str(removal_time) + str(cask_name)),
            masses["I129"],
            atomic_masses["I129"],
            half_lives["I129"],
            removal_time,
            149,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Np239"],
            ("Np239" + str(removal_time) + str(cask_name)),
            masses["Np239"],
            atomic_masses["Np239"],
            half_lives["Np239"],
            removal_time,
            714,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Tc99"],
            ("Tc99" + str(removal_time) + str(cask_name)),
            masses["Tc99"],
            atomic_masses["Tc99"],
            half_lives["Tc99"],
            removal_time,
            440,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Zr93"],
            ("Zr93" + str(removal_time) + str(cask_name)),
            masses["Zr93"],
            atomic_masses["Zr93"],
            half_lives["Zr93"],
            removal_time,
            1000,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Ce144"],
            ("Ce144" + str(removal_time) + str(cask_name)),
            masses["Ce144"],
            atomic_masses["Ce144"],
            half_lives["Ce144"],
            removal_time,
            318,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Kr88"],
            ("Kr88" + str(removal_time) + str(cask_name)),
            masses["Kr88"],
            atomic_masses["Kr88"],
            half_lives["Kr88"],
            removal_time,
            2918,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Pr144"],
            ("Pr144" + str(removal_time) + str(cask_name)),
            masses["Pr144"],
            atomic_masses["Pr144"],
            half_lives["Pr144"],
            removal_time,
            2997,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Rb88"],
            ("Rb88" + str(removal_time) + str(cask_name)),
            masses["Rb88"],
            atomic_masses["Rb88"],
            half_lives["Rb88"],
            removal_time,
            5301,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Rh106"],
            ("Rh106" + str(removal_time) + str(cask_name)),
            masses["Rh106"],
            atomic_masses["Rh106"],
            half_lives["Rh106"],
            removal_time,
            3541,
        ),
    )
    spectra.Add(
        load_spec(
            isotope_data["Ru106"],
            ("Ru106" + str(removal_time) + str(cask_name)),
            masses["Ru106"],
            atomic_masses["Ru106"],
            half_lives["Ru106"],
            removal_time,
            39,
        ),
    )

    # adding on the spectrum of any newly created isotopes from previous decays
    # all of these decay chains have a branching ratio of 1
    # if any additional isotopes were to be added with decay chains involving more beta emitting isotopes then they can be added here

    if removal_time != 0:
        extra_Y90 = AddDecays(
            t=removal_time,
            parent_prop0=Sr90_prop,
            parent_half_life=half_lives["Sr90"],
            daughter_half_life=half_lives["Y90"],
            total_mass=total_m,
        )
        spectra.Add(
            load_spec(
                isotope_data["Y90"],
                "additional Y90" + str(removal_time) + str(cask_name),
                extra_Y90,
                atomic_masses["Y90"],
                half_lives["Y90"],
                0,
                2278,
            ),
        )

        extra_Pr144 = AddDecays(
            t=removal_time,
            parent_prop0=Ce144_prop,
            parent_half_life=half_lives["Ce144"],
            daughter_half_life=half_lives["Pr144"],
            total_mass=total_m,
        )
        spectra.Add(
            load_spec(
                isotope_data["Pr144"],
                "additional Pr144" + str(removal_time) + str(cask_name),
                extra_Pr144,
                atomic_masses["Pr144"],
                half_lives["Pr144"],
                0,
                2997,
            ),
        )

        extra_Rb88 = AddDecays(
            t=removal_time,
            parent_prop0=Kr88_prop,
            parent_half_life=half_lives["Kr88"],
            daughter_half_life=half_lives["Rb88"],
            total_mass=total_m,
        )
        spectra.Add(
            load_spec(
                isotope_data["Rb88"],
                "additional Rb88" + str(removal_time) + str(cask_name),
                extra_Rb88,
                atomic_masses["Rb88"],
                half_lives["Rb88"],
                0,
                5301,
            ),
        )

        extra_Rh106 = AddDecays(
            t=removal_time,
            parent_prop0=Ru106_prop,
            parent_half_life=half_lives["Ru106"],
            daughter_half_life=half_lives["Rh106"],
            total_mass=total_m,
        )
        spectra.Add(
            load_spec(
                isotope_data["Rh106"],
                "additional Rh106" + str(removal_time) + str(cask_name),
                extra_Rh106,
                atomic_masses["Rh106"],
                half_lives["Rh106"],
                0,
                3541,
            ),
        )

    total_spec = add_spec(spectra)
    total_spec.SetTitle("Total Spectrum")
    total_spec.GetXaxis().SetRangeUser(0, max_E)
    return total_spec
