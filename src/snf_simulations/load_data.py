from importlib.resources import files

import numpy as np


def load_antineutrino_data():
    """Load in antineutrino spectrum data from IAEA."""
    spec_files = files("snf_simulations.data.spec_data")

    # Columns represent energy, dN/dE and e=uncertainty on flux respectively.
    # Some elements like Y90 have two decay chains, so they are cut off where the
    # main decay chain ends.

    Sr90 = np.genfromtxt(spec_files.joinpath("Sr90_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]

    Y90 = np.genfromtxt(spec_files.joinpath("Y90_an.txt"), skip_header=1)
    end = int(np.where(Y90[:, 7] == 2278.5)[0])
    Y90 = Y90[:end, [7, 10, 11]]

    Pu241 = np.genfromtxt(spec_files.joinpath("Pu241_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]
    Cs137 = np.genfromtxt(spec_files.joinpath("Cs137_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]
    Am242 = np.genfromtxt(spec_files.joinpath("Am242_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]
    Cs135 = np.genfromtxt(spec_files.joinpath("Cs135_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]
    I129 = np.genfromtxt(spec_files.joinpath("I129_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]
    Np239 = np.genfromtxt(spec_files.joinpath("Np239_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]

    Tc99 = np.genfromtxt(spec_files.joinpath("Tc99_an.txt"), skip_header=1)
    end0 = int(np.where(Tc99[:, 7] == 297.5)[0])
    Tc99 = Tc99[:end0, [7, 10, 11]]

    Zr93 = np.genfromtxt(spec_files.joinpath("Zr93_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]

    Ce144 = np.genfromtxt(spec_files.joinpath("Ce144_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]
    Kr88 = np.genfromtxt(spec_files.joinpath("Kr88_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]
    Pr144 = np.genfromtxt(spec_files.joinpath("Pr144_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]
    Rb88 = np.genfromtxt(spec_files.joinpath("Rb88_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]
    Rh106 = np.genfromtxt(spec_files.joinpath("Rh106_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]
    Ru106 = np.genfromtxt(spec_files.joinpath("Ru106_an.txt"), skip_header=1)[
        :, [7, 10, 11],
    ]

    return (
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
    )
