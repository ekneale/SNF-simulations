"""Module for loading FISPIN .tbQ output files."""

from pathlib import Path

import numpy as np
import pandas as pd

from .utils import _UNITS_TO_SECONDS

_UNITS_TO_YEARS = {
    "SECONDS": _UNITS_TO_SECONDS["sec"] / _UNITS_TO_SECONDS["year"],
    "MINS": _UNITS_TO_SECONDS["minute"] / _UNITS_TO_SECONDS["year"],
    "HOURS": _UNITS_TO_SECONDS["hour"] / _UNITS_TO_SECONDS["year"],
    "DAYS": _UNITS_TO_SECONDS["day"] / _UNITS_TO_SECONDS["year"],
    "YEARS": 1.0,
}


def load_tabqfile(filepath_or_contents: str | Path) -> dict[str, pd.DataFrame]:
    """Load in a FISPIN .tbQ output file and extract the data.

    Args:
        filepath_or_contents: Path to the file to load,
        or the contents of the file as a string.

    Returns:
        Dictionary of pandas dataframes containing the isotope data,
        for each time step in the file.
        Keys are the simulation time string values from the file,
        e.g. "6.000E+01 MINS" or "2.800E+01 DAYS".

    """
    try:
        with open(filepath_or_contents) as f:
            lines = f.readlines()
    except FileNotFoundError:
        if isinstance(filepath_or_contents, Path):
            raise
        # If the input is not a valid file path, treat it as the file contents.
        lines = filepath_or_contents.splitlines(keepends=True)

    # The file can contain multiple sections for different time steps.
    # Header format is "*** TIME    2.347E+00 YEARS"
    # We'll extract the value and unit string, but leave converting it to the user.
    time_dfs = {}
    for i, line in enumerate(lines):
        if not line.startswith("*** TIME"):
            continue

        # Get the time value and unit from the header line
        time_str = " ".join(line.split()[2:4])

        # The next line should be the header
        header_keys = lines[i + 1].split()

        # Read the following lines until we hit the total at the end of the section
        # (or the end of the file)
        isotope_data = []
        for data_line in lines[i + 2 :]:
            if data_line.startswith("TOTAL"):
                break
            # Complication: the first few characters are the isotope,
            # which can contain spaces (e.g. "U 235").
            # We'll remove spaces from the first 12 characters before splitting.
            # We also format the isotope to be capitalised (e.g. "Sr90"),
            # by default the file has them in uppercase (e.g. "SR90").
            isotope = data_line[0:12].capitalize().replace(" ", "")
            data_line_new = f"{isotope:>12}{data_line[12:]}"
            data = data_line_new.strip().split()
            isotope_data.append(data)

        # Convert to a DataFrame and store
        df = pd.DataFrame(isotope_data, columns=pd.Index(header_keys))
        for col in header_keys[1:]:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        time_dfs[time_str] = df

    return time_dfs


def _convert_sim_time_to_years(time_str: str) -> float:
    """Convert a simulation time string from the .tbQ file into years."""
    value, unit = time_str.split()
    return float(value) * _UNITS_TO_YEARS[unit]


def get_isotope_masses(
    filepath: str | Path, time_str: str | None = None
) -> tuple[dict[str, float], float]:
    """Get the isotope masses from a FISPIN .tbQ output file.

    Args:
        filepath: Path to the file to load.
        time_str: Specific simulation time to extract data for.
            If None, uses the earliest time in the file.
            Will raise an error if the specified string is not found in the file.

    Returns:
        isotope_masses: Dictionary of isotope masses at the specified
            time, where keys are isotope names and values are the mass of each
            isotope in kg, e.g. {"Sr90": 5.356e-4, "Y90": 1.3922e-7, ...}.
        cooling_time: The cooling time in years that the masses correspond to,
            converted from the file time units.

    """
    time_dfs = load_tabqfile(filepath)

    if time_str is not None:
        # Select the dataframe for the specified time, and convert it to years.
        if time_str not in time_dfs:
            msg = f"Specified time string '{time_str}' not found in file: {filepath}"
            msg += f"\nAvailable times: {list(time_dfs.keys())}"
            raise ValueError(msg)
        df = time_dfs[time_str]
        cooling_time = _convert_sim_time_to_years(time_str)
    else:
        # If no time string is specified, take the isotope masses from the earliest
        # simulated cooling time.
        selected_time_str = ""
        cooling_time = np.inf
        for inner_time_str in time_dfs:
            time_years = _convert_sim_time_to_years(inner_time_str)
            if time_years < cooling_time:
                cooling_time = time_years
                selected_time_str = inner_time_str
        df = time_dfs[selected_time_str]

    # Return the masses of each isotope and the selected cooling time.
    masses = pd.to_numeric(df["GRAMS"], errors="coerce")
    names = df["ALL-NUC"].astype(str)
    isotope_masses = {
        name: mass * 1e-3  # convert from grams to kg
        for name, mass in zip(names, masses, strict=True)
    }
    return isotope_masses, cooling_time
