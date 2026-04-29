"""Module for loading reactor data."""

from pathlib import Path

import pandas as pd

_UNITS_TO_SECONDS = {
    "SECONDS": 1,
    "MINS": 60,
    "HOURS": 60 * 60,
    "DAYS": 24 * 60 * 60,
    "YEARS": 365 * 24 * 60 * 60,
}


def load_tabqfile(filepath: str | Path) -> dict[str, pd.DataFrame]:
    """Load in a FISPIN .tbQ file and extract the data.

    Args:
        filepath: Path to the .tabq file to load.

    Returns:
        Dictionary of pandas dataframes containing the isotope data,
        for each time step in the file.
        Keys are the time values from the file,
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


def get_isotope_proportions(
    filepath: str | Path, timestep: str | None = None
) -> dict[str, float]:
    """Get the isotope proportions from a FISPIN .tabq file.

    Args:
        filepath: Path to the .tabq file to load.
        timestep: Specific time step to extract data for.
            If None, uses the earliest time step in the file.

    Returns:
        Dictionary of isotope proportions in the fuel,
        e.g. {"Sr90": 5.356e-4, "Y90": 1.3922e-7, ...}.

    """
    time_dfs = load_tabqfile(filepath)
    if timestep is None:
        # Return the isotope proportions from the earliest time step
        # Annoyingly, the file doesn't have to be in order and the keys are strings
        # in different units, so we have to convert them here to find the earliest.
        try:
            timesteps = {
                (
                    float(time_str.split()[0]) * _UNITS_TO_SECONDS[time_str.split()[1]]
                ): time_str
                for time_str in time_dfs
            }
        except (IndexError, KeyError) as err:
            msg = f"Error parsing time steps from file: {filepath}"
            raise ValueError(msg) from err
        timestep = timesteps[min(timesteps)]

    # Calculate the proportions of each isotope of the total
    df = time_dfs[timestep]
    total_mass = df["GRAMS"].sum()
    isotope_masses = df.set_index("ALL-NUC")["GRAMS"].to_dict()
    proportions = {
        isotope: float(mass / total_mass) for isotope, mass in isotope_masses.items()
    }
    return proportions
