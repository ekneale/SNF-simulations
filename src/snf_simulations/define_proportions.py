import numpy as np
import ROOT

from .load_data import load_antineutrino_data, load_isotopes
from .spec import add_spec, load_spec


def AddDecays(
    t, parent_prop0, parent_half_life, daughter_half_life, total_mass=1000, BR=1,
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

    # Load the isotope data
    isotopes, atomic_masses, half_lives = load_isotopes()

    (
        Sr90,
        Y90,
        Pu241,
        Cs137,
        Am242,
        Cs135,
        I129,
        Np239,
        Tc99,
        Zr93,
        Ce144,
        Kr88,
        Pr144,
        Rb88,
        Rh106,
        Ru106,
    ) = load_antineutrino_data(isotopes)


    # inputted proportions of isotopes defined by the user
    proportions = [
        Sr90_prop,
        Y90_prop,
        Pu241_prop,
        Cs137_prop,
        Am242_prop,
        Cs135_prop,
        I129_prop,
        Np239_prop,
        Tc99_prop,
        Zr93_prop,
        Ce144_prop,
        Kr88_prop,
        Pr144_prop,
        Rb88_prop,
        Rh106_prop,
        Ru106_prop,
    ]

    # calculated mass of isotopes from inputted proportions
    masses = [prop * total_m for prop in proportions]

    spectra = ROOT.TList()

    # adding the scaled spectra of each isotope to the list that will later be added together

    spectra.Add(
        load_spec(
            Sr90,
            (isotopes[0] + str(removal_time) + str(cask_name)),
            masses[0],
            atomic_masses[0],
            half_lives[0],
            removal_time,
            546,
        ),
    )
    spectra.Add(
        load_spec(
            Y90,
            (isotopes[1] + str(removal_time) + str(cask_name)),
            masses[1],
            atomic_masses[1],
            half_lives[1],
            removal_time,
            2278,
        ),
    )
    spectra.Add(
        load_spec(
            Pu241,
            (isotopes[2] + str(removal_time) + str(cask_name)),
            masses[2],
            atomic_masses[2],
            half_lives[2],
            removal_time,
            20,
        ),
    )
    spectra.Add(
        load_spec(
            Cs137,
            (isotopes[3] + str(removal_time) + str(cask_name)),
            masses[3],
            atomic_masses[3],
            half_lives[3],
            removal_time,
            1175,
        ),
    )
    spectra.Add(
        load_spec(
            Am242,
            (isotopes[4] + str(removal_time) + str(cask_name)),
            masses[4],
            atomic_masses[4],
            half_lives[4],
            removal_time,
            664,
        ),
    )
    spectra.Add(
        load_spec(
            Cs135,
            (isotopes[5] + str(removal_time) + str(cask_name)),
            masses[5],
            atomic_masses[5],
            half_lives[5],
            removal_time,
            268,
        ),
    )
    spectra.Add(
        load_spec(
            I129,
            (isotopes[6] + str(removal_time) + str(cask_name)),
            masses[6],
            atomic_masses[6],
            half_lives[6],
            removal_time,
            149,
        ),
    )
    spectra.Add(
        load_spec(
            Np239,
            (isotopes[7] + str(removal_time) + str(cask_name)),
            masses[7],
            atomic_masses[7],
            half_lives[7],
            removal_time,
            714,
        ),
    )
    spectra.Add(
        load_spec(
            Tc99,
            (isotopes[8] + str(removal_time) + str(cask_name)),
            masses[8],
            atomic_masses[8],
            half_lives[8],
            removal_time,
            440,
        ),
    )
    spectra.Add(
        load_spec(
            Zr93,
            (isotopes[9] + str(removal_time) + str(cask_name)),
            masses[9],
            atomic_masses[9],
            half_lives[9],
            removal_time,
            1000,
        ),
    )
    spectra.Add(
        load_spec(
            Ce144,
            (isotopes[10] + str(removal_time) + str(cask_name)),
            masses[10],
            atomic_masses[10],
            half_lives[10],
            removal_time,
            318,
        ),
    )
    spectra.Add(
        load_spec(
            Kr88,
            (isotopes[11] + str(removal_time) + str(cask_name)),
            masses[11],
            atomic_masses[11],
            half_lives[11],
            removal_time,
            2918,
        ),
    )
    spectra.Add(
        load_spec(
            Pr144,
            (isotopes[12] + str(removal_time) + str(cask_name)),
            masses[12],
            atomic_masses[12],
            half_lives[12],
            removal_time,
            2997,
        ),
    )
    spectra.Add(
        load_spec(
            Rb88,
            (isotopes[13] + str(removal_time) + str(cask_name)),
            masses[13],
            atomic_masses[13],
            half_lives[13],
            removal_time,
            5301,
        ),
    )
    spectra.Add(
        load_spec(
            Rh106,
            (isotopes[14] + str(removal_time) + str(cask_name)),
            masses[14],
            atomic_masses[14],
            half_lives[14],
            removal_time,
            3541,
        ),
    )
    spectra.Add(
        load_spec(
            Ru106,
            (isotopes[15] + str(removal_time) + str(cask_name)),
            masses[15],
            atomic_masses[15],
            half_lives[15],
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
            parent_half_life=half_lives[0],
            daughter_half_life=half_lives[1],
            total_mass=total_m,
        )
        spectra.Add(
            load_spec(
                Y90,
                "additional Y90" + str(removal_time) + str(cask_name),
                extra_Y90,
                atomic_masses[1],
                half_lives[1],
                0,
                2278,
            ),
        )

        extra_Pr144 = AddDecays(
            t=removal_time,
            parent_prop0=Ce144_prop,
            parent_half_life=half_lives[10],
            daughter_half_life=half_lives[12],
            total_mass=total_m,
        )
        spectra.Add(
            load_spec(
                Pr144,
                "additional Pr144" + str(removal_time) + str(cask_name),
                extra_Pr144,
                atomic_masses[12],
                half_lives[12],
                0,
                2997,
            ),
        )

        extra_Rb88 = AddDecays(
            t=removal_time,
            parent_prop0=Kr88_prop,
            parent_half_life=half_lives[11],
            daughter_half_life=half_lives[13],
            total_mass=total_m,
        )
        spectra.Add(
            load_spec(
                Rb88,
                "additional Rb88" + str(removal_time) + str(cask_name),
                extra_Rb88,
                atomic_masses[13],
                half_lives[13],
                0,
                5301,
            ),
        )

        extra_Rh106 = AddDecays(
            t=removal_time,
            parent_prop0=Ru106_prop,
            parent_half_life=half_lives[15],
            daughter_half_life=half_lives[14],
            total_mass=total_m,
        )
        spectra.Add(
            load_spec(
                Rh106,
                "additional Rh106" + str(removal_time) + str(cask_name),
                extra_Rh106,
                atomic_masses[14],
                half_lives[14],
                0,
                3541,
            ),
        )

    total_spec = add_spec(spectra)
    total_spec.SetTitle("Total Spectrum")
    total_spec.GetXaxis().SetRangeUser(0, max_E)
    return total_spec
