"""Load in data for antineutrino spectra from IAEA files."""

from importlib.resources import files

import numpy as np


def load_isotope_data(isotopes=None):
    """Load in isotope parameters from the database.

    Args:
        isotopes (list of str, optional): List of isotope names to load data for.
            If None, load data for all isotopes in the database.

    Returns:
        tuple: Two dictionaries:
            - atomic_masses: Dictionary of atomic masses (g/mol) for each isotope.
            - half_lives: Dictionary of half-lives (years) for each isotope.

    """
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
    if isotopes is not None:
        data = data[np.isin(data[:, 0], isotopes)]
    atomic_masses = {str(d[0]): int(d[1]) for d in data}
    half_lives = {str(d[0]): float(d[2]) for d in data}
    return atomic_masses, half_lives


def load_spectrum(isotope_name):
    """Load in antineutrino spectrum data from text files.

    Args:
        isotope_name (str): Name of the isotope to load data for.

    Returns:
        np.ndarray: Array containing energy, dN/dE, and uncertainty.

    """
    # TODO: download spectra from IAEA database, and cache locally
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
        dict of np.ndarray: Dictionary of arrays containing spectrum data
        for each isotope.
        Data arrays contain energy, dN/dE, and uncertainty.

    """
    data = {}
    for isotope in isotopes:
        isotope_data = load_spectrum(isotope)
        data[isotope] = isotope_data
    return data
