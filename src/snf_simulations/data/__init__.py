"""Data loading module for antineutrino spectra calculations.

Provides functions to load isotope data, reactor fuel composition data, and
antineutrino spectrum data from the built-in database files.
"""

from importlib.resources import as_file, files

import numpy as np


def get_reactors() -> list[str]:
    """Get a list of available reactors from the database.

    Returns:
        List of reactor names.

    """
    data_files = files("snf_simulations.data.reactor_data")
    with as_file(data_files) as path:
        return [file.stem for file in path.iterdir() if file.suffix == ".csv"]


def load_reactor_data(reactor: str) -> dict[str, float]:
    """Load in reactor isotope proportions from the database.

    Args:
        reactor: Name of the reactor to load data for.

    Returns:
        Dictionary of isotope proportions in the fuel,
        e.g. {"Sr90": 5.356e-4, "Y90": 1.3922e-7, ...}.

    """
    data_files = files("snf_simulations.data.reactor_data")
    filename = data_files / f"{reactor}.csv"

    # Check if file exists
    with as_file(filename) as filepath:
        if not filepath.is_file():
            msg = f"Reactor {reactor} data file not found."
            valid_reactors = ", ".join(get_reactors())
            msg += f" Valid reactors are: {valid_reactors}."
            raise ValueError(msg)

        data = np.genfromtxt(
            filepath,
            delimiter=",",
            skip_header=1,
            dtype=str,
        )

    return {str(d[0]): float(d[1]) for d in data}


def load_isotope_data(
    isotopes: list[str] | None = None,
) -> tuple[dict[str, int], dict[str, float]]:
    """Load in isotope parameters from the database.

    Args:
        isotopes: List of isotope names to load data for.
            If None, load data for all isotopes in the database.

    Returns:
        Tuple of two dictionaries:
        - molar_masses: Dictionary of molar masses (g/mol) for each isotope.
        - half_lives: Dictionary of half-lives (years) for each isotope.

    """
    data_files = files("snf_simulations.data")
    filename = data_files / "isotopes.csv"

    with as_file(filename) as filepath:
        if not filepath.is_file():
            msg = "Isotope CSV file not found."
            raise ValueError(msg)

        data = np.genfromtxt(
            filepath,
            delimiter=",",
            skip_header=1,
            dtype=str,
        )

    if isotopes is not None:
        data = data[np.isin(data[:, 0], isotopes)]
    molar_masses = {str(d[0]): int(d[1]) for d in data}
    half_lives = {str(d[0]): float(d[2]) for d in data}
    return molar_masses, half_lives


def load_spectrum(isotope_name: str) -> np.ndarray:
    """Load in antineutrino spectrum data from text files.

    Args:
        isotope_name: Name of the isotope to load data for.

    Returns:
        Array containing energy, dN/dE, and uncertainty.

    """
    # TODO: download spectra from IAEA database, and cache locally
    spec_files = files("snf_simulations.data.spec_data")
    filename = spec_files / f"{isotope_name}_an.txt"

    with as_file(filename) as filepath:
        if not filepath.is_file():
            msg = f"Spectrum data file for {isotope_name} not found."
            raise ValueError(msg)

        data = np.genfromtxt(filepath, skip_header=1)

    # Some isotopes have multiple decay chains, so cut off where the
    # main decay chain ends based on the p_energy column.
    data = data[data[:, 3] == 0]

    return data[:, [7, 10, 11]]  # energy, dN/dE, uncertainty


def load_antineutrino_data(isotopes: list[str]) -> dict[str, np.ndarray]:
    """Load in antineutrino spectrum data from the IAEA.

    Args:
        isotopes: List of isotope names to load data for.

    Returns:
        Dictionary of arrays containing spectrum data for each isotope.
        Data arrays contain energy, dN/dE, and uncertainty.

    """
    data = {}
    for isotope in isotopes:
        isotope_data = load_spectrum(isotope)
        data[isotope] = isotope_data
    return data
