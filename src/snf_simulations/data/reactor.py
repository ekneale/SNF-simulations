"""Module for loading reactor data."""

from importlib.resources import as_file, files

import numpy as np


def get_reactors() -> list[str]:
    """Get a list of available reactors from the database.

    Returns:
        List of reactor names.

    """
    data_files = files("snf_simulations.data") / "reactor_data"
    with as_file(data_files) as path:
        return [file.stem for file in path.iterdir() if file.suffix == ".csv"]


def get_reactor_data(reactor: str) -> dict[str, float]:
    """Load in reactor isotope proportions from the database.

    Args:
        reactor: Name of the reactor to load data for.

    Returns:
        Dictionary of isotope proportions in the fuel,
        e.g. {"Sr90": 5.356e-4, "Y90": 1.3922e-7, ...}.

    """
    filename = files("snf_simulations.data") / "reactor_data" / f"{reactor}.csv"

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
