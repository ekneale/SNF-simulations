"""Data loading module for antineutrino spectra calculations.

Provides functions to load isotope data, reactor fuel composition data, and
antineutrino spectrum data from the built-in database files.
"""

import os
import re
import urllib.request
from importlib.resources import as_file, files
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

_CACHE_DIR_ENV_VAR = "SNF_SIMULATIONS_CACHE_DIR"


def get_reactors() -> list[str]:
    """Get a list of available reactors from the database.

    Returns:
        List of reactor names.

    """
    data_files = files("snf_simulations.data") / "reactor_data"
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


def _parse_nuclide(nuclide: str) -> str:
    """Parse a nuclide name into the format expected by the IAEA database.

    Args:
        nuclide: Nuclide name in the format 'ElementMass' or 'MassElement',
            e.g. 'Ru106' or '106Ru'.

    Returns:
        Nuclide name in the format 'masselement', e.g. '106ru'.
        Note that the element symbol is converted to lowercase.

    """
    match = re.match(r"^([A-Za-z]+)(\d+)$", nuclide)
    if match:
        element, mass = match.groups()
        return f"{mass}{element.lower()}"
    match = re.match(r"^(\d+)([A-Za-z]+)$", nuclide)
    if match:
        mass, element = match.groups()
        return f"{mass}{element.lower()}"
    else:
        msg = f"Nuclide format not recognized: {nuclide}"
        msg += " Use 'ElementMass' or 'MassElement', e.g., 'Ru106' or '106Ru'."
        raise ValueError(msg)


def _get_cache_dir() -> Path:
    """Return the writable directory used for downloaded spectrum data."""
    cache_dir = os.environ.get(_CACHE_DIR_ENV_VAR)
    if cache_dir is not None:
        path = Path(cache_dir).expanduser()
    else:
        xdg_cache_home = os.environ.get("XDG_CACHE_HOME")
        if xdg_cache_home is not None:
            path = Path(xdg_cache_home).expanduser() / "snf_simulations" / "spec_data"
        else:
            path = Path.home() / ".cache" / "snf_simulations" / "spec_data"

    path.mkdir(parents=True, exist_ok=True)
    return path


def _get_cache_file(nuclide: str) -> Path:
    """Return the writable cache file path for a nuclide."""
    return _get_cache_dir() / f"{nuclide}.csv"


def download_spectrum_data(nuclide: str) -> str:
    """Download the antineutrino spectrum for a given nuclide from the IAEA database.

    Args:
        nuclide: Name of the nuclide to get the spectrum for, e.g. "Ru106" or "106Ru".

    Returns:
        Path to the downloaded spectrum data file in the cache.

    """
    # Download the data file
    nuclide = _parse_nuclide(nuclide)
    url = (
        "https://nds.iaea.org/relnsd/v1/data?fields=bin_beta"
        f"&nuclides={nuclide}&rad_types=bm"
    )
    req = urllib.request.Request(url)  # noqa: S310
    req.add_header(
        "User-Agent",
        "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox77.0",
    )
    try:
        content = urllib.request.urlopen(req, timeout=5).read().decode("utf-8")  # noqa: S310
    except Exception as err:
        msg = f"Error downloading spectrum data for {nuclide} from IAEA database: {err}"
        raise RuntimeError(msg) from err

    data = pd.read_csv(StringIO(content))

    # Some nuclides (e.g. Ru106) have duplicate rows in the IAEA database.
    # These break creating histograms, so remove exact duplicates before caching.
    data = data.drop_duplicates()

    filename = _get_cache_file(nuclide)
    if not filename.is_file():
        data.to_csv(filename, index=False)

    return str(filename)


def load_spectrum_file(nuclide: str) -> np.ndarray:
    """Load in antineutrino spectrum data from a CSV file in the cache.

    Args:
        nuclide: Name of the nuclide to load data for.

    Returns:
        Array containing energy, flux, and uncertainty.

    """
    nuclide = _parse_nuclide(nuclide)
    cache_file = _get_cache_file(nuclide)
    if not cache_file.is_file():
        msg = f"Spectrum data file for {nuclide} not found in cache."
        raise ValueError(msg)

    df = pd.read_csv(cache_file)

    # Some isotopes have multiple decay chains, so cut off where the
    # main decay chain ends based on the p_energy column.
    df = df[df["p_energy"] == 0]

    # Return only the relevant columns as a "D numpy array.
    return df[["bin_en", "dn_de_nu", "unc_dn_de_nu"]].to_numpy()


def load_spectrum(nuclide: str) -> np.ndarray:
    """Load in antineutrino spectrum data for a given nuclide.

    If the spectrum data is not already in the cache, it is downloaded from the
    IAEA database and saved locally before loading.

    Args:
        nuclide: Name of the nuclide to load data for.

    Returns:
        Array containing energy, flux, and uncertainty.

    """
    nuclide = _parse_nuclide(nuclide)
    if not _get_cache_file(nuclide).is_file():
        download_spectrum_data(nuclide)
    return load_spectrum_file(nuclide)


def load_antineutrino_data(nuclides: list[str]) -> dict[str, np.ndarray]:
    """Load in IAEA antineutrino spectrum data for the specified nuclides.

    Args:
        nuclides: List of nuclide names to load data for.

    Returns:
        Dictionary of arrays containing spectrum data for each nuclide.
        Data arrays contain energy, flux, and uncertainty.

    """
    data = {}
    for nuclide in nuclides:
        nuclide_data = load_spectrum(nuclide)
        data[nuclide] = nuclide_data
    return data
