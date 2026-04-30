"""Module for loading reactor data."""

from pathlib import Path

import numpy as np
import pandas as pd

_UNITS_TO_SECONDS = {
    "SECONDS": 1,
    "MINS": 60,
    "HOURS": 60 * 60,
    "DAYS": 24 * 60 * 60,
    "YEARS": 365 * 24 * 60 * 60,
}


def load_tabqfile(filepath: str | Path) -> dict[str, pd.DataFrame]:
    """Load in a FISPIN .tbQ output file and extract the data.

    Args:
        filepath: Path to the file to load.

    Returns:
        Dictionary of pandas dataframes containing the isotope data,
        for each time step in the file.
        Keys are the simulation time string values from the file,
        e.g. "6.000E+01 MINS" or "2.800E+01 DAYS".

    """
    with open(filepath) as f:
        lines = f.readlines()

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
    """Convert a simulation time string from the .tabq file into years."""
    value, unit = time_str.split()
    return float(value) * _UNITS_TO_SECONDS[unit] / _UNITS_TO_SECONDS["YEARS"]


def get_isotope_proportions(
    filepath: str | Path, time_str: str | None = None
) -> tuple[dict[str, float], float]:
    """Get the isotope proportions from a FISPIN .tabq file.

    Args:
        filepath: Path to the .tabq file to load.
        time_str: Specific simulation time to extract data for.
            If None, uses the earliest time in the file.
            Will raise an error if the specified string is not found in the file.

    Returns:
        isotope_proportions: Dictionary of isotope proportions at the specified
            time, where keys are isotope names and values are the proportion of the
            total mass for each isotope, e.g. {"Sr90": 5.356e-4, "Y90": 1.3922e-7, ...}.
        cooling_time: The cooling time in years that the proportions correspond to,
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
        # If no time string is specified, use the isotope proportions from the earliest
        # simulated cooling time.
        selected_time_str = ""
        cooling_time = np.inf
        for inner_time_str in time_dfs:
            time_years = _convert_sim_time_to_years(inner_time_str)
            if time_years < cooling_time:
                cooling_time = time_years
                selected_time_str = inner_time_str
        df = time_dfs[selected_time_str]

    # Calculate the proportions of each isotope of the total mass.
    isotope_masses = pd.to_numeric(df["GRAMS"], errors="coerce")
    total_mass = float(isotope_masses.sum())
    isotope_names = df["ALL-NUC"].astype(str)
    proportions: dict[str, float] = {}
    for isotope, mass in zip(isotope_names, isotope_masses, strict=True):
        proportions[str(isotope)] = float(mass) / total_mass

    return proportions, cooling_time
