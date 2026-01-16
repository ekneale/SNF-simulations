"""Load in data for antineutrino spectra from IAEA files."""
from importlib.resources import files

import numpy as np


def load_isotopes():
    """Load in isotope parameters from the database."""
    data_files = files("snf_simulations.data")
    filename = data_files.joinpath("isotopes.csv")
    if not filename.is_file():
        msg = "Isotope CSV file not found."
        raise ValueError(msg)

    data = np.genfromtxt(
        filename,
        delimiter=",",
        skip_header=1,
        dtype=str,
    )
    names, atomic_masses, half_lives = data[:,0], data[:,1], data[:,2]
    atomic_masses = atomic_masses.astype(float)
    half_lives = half_lives.astype(float)
    return names, atomic_masses, half_lives


def load_spec_data(isotope_name):
    """Load in spectrum data from text files.

    Args:
        isotope_name (str): Name of the isotope to load data for.

    Returns:
        np.ndarray: Array containing energy, dN/dE, and uncertainty.

    """
    spec_files = files("snf_simulations.data.spec_data")
    filename = spec_files.joinpath(f"{isotope_name}_an.txt")
    if not filename.is_file():
        msg = f"Spectrum data file for {isotope_name} not found."
        raise ValueError(msg)

    data = np.genfromtxt(spec_files.joinpath(filename), skip_header=1)

    # Some isotopes have multiple decay chains, so cut off where the
    # main decay chain ends based on the p_energy column.
    data = data[data[:, 3] == 0]

    return data[:, [7, 10, 11]]  # energy, dN/dE, uncertainty


def load_antineutrino_data(isotopes):
    """Load in antineutrino spectrum data from the IAEA.

    Args:
        isotopes (list of str): List of isotope names to load data for.

    Returns:
        list of np.ndarray: List of arrays containing spectrum data for each isotope.

    """
    spec_data = []
    for isotope in isotopes:
        data = load_spec_data(isotope)
        spec_data.append(data)
    return spec_data
