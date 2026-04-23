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
from mendeleev import isotope

_CACHE_DIR_ENV_VAR = "SNF_SIMULATIONS_CACHE_DIR"


# mendeleev uses this conversion, see https://github.com/lmmentel/mendeleev/pull/160
_SECONDS_PER_YEAR = 31_556_926.0
_UNITS_TO_SECONDS = {
    "ysec": 1e-24,
    "zsec": 1e-21,
    "asec": 1e-18,
    "psec": 1e-12,
    "nsec": 1e-9,
    "usec": 1e-6,
    "msec": 1e-3,
    "sec": 1.0,
    "minute": 60.0,
    "hour": 3600.0,
    "day": 86400.0,
    "year": _SECONDS_PER_YEAR,
    "kyear": 1e3 * _SECONDS_PER_YEAR,
    "Myear": 1e6 * _SECONDS_PER_YEAR,
    "Gyear": 1e9 * _SECONDS_PER_YEAR,
    "Tyear": 1e12 * _SECONDS_PER_YEAR,
    "Pyear": 1e15 * _SECONDS_PER_YEAR,
    "Eyear": 1e18 * _SECONDS_PER_YEAR,
    "Zyear": 1e21 * _SECONDS_PER_YEAR,
    "Yyear": 1e24 * _SECONDS_PER_YEAR,
}


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


def _parse_isotope(isotope_name: str) -> tuple[str, int]:
    """Parse an isotope name and return its element and mass number.

    Args:
        isotope_name: Name of the isotope to parse.
            Format should be 'ElementMass' (e.g. Ru106) or 'MassElement' (e.g. 106Ru).

    Returns:
        element: The element symbol (case unchanged).
        mass_number: The mass number.

    """
    match = re.match(r"^([A-Za-z]+)(\d+)$", isotope_name)
    if match:
        element, mass_number = match.groups()
        return element, int(mass_number)
    match = re.match(r"^(\d+)([A-Za-z]+)$", isotope_name)
    if match:
        mass_number, element = match.groups()
        return element, int(mass_number)
    else:
        msg = f"Isotope format not recognized: {isotope_name}"
        msg += " Use 'ElementMass' or 'MassElement', e.g., 'Ru106' or '106Ru'."
        raise ValueError(msg)


def get_isotope_properties(isotope_name: str) -> dict[str, float]:
    """Get the mass and half-life for the given isotope using the mendeleev package.

    Args:
        isotope_name: Name of the isotope.
            Format should be 'ElementMass' (e.g. Ru106) or 'MassElement' (e.g. 106Ru).

    Returns:
        Dictionary containing the molar mass (in g/mol) and half-life (in years)
        of the isotope.

    """
    element, mass_number = _parse_isotope(isotope_name)
    mendeleev_isotope = isotope(element, mass_number)

    molar_mass = float(mendeleev_isotope.mass)  # ty: ignore
    half_life = float(mendeleev_isotope.half_life)  # ty: ignore
    unit = str(mendeleev_isotope.half_life_unit)
    seconds_per_unit = _UNITS_TO_SECONDS.get(unit)
    if seconds_per_unit is None:
        msg = f"Unsupported half-life unit for isotope {isotope_name}: {unit}. "
        msg += "Supported units are: "
        msg += ", ".join(sorted(_UNITS_TO_SECONDS))
        raise ValueError(msg)
    half_life_years = half_life * seconds_per_unit / _SECONDS_PER_YEAR

    return {"molar_mass": molar_mass, "half_life": half_life_years}


def _get_cache_dir() -> Path:
    """Return the writable directory used for downloaded spectrum data.

    Returns:
        Path to the cache directory.

    """
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


def _get_cache_file(isotope_name: str) -> Path:
    """Return the writable cache file path for an isotope.

    Args:
        isotope_name: Name of the isotope to get the cache file for.
            Format should be 'ElementMass' (e.g. Ru106) or 'MassElement' (e.g. 106Ru).

    Returns:
        Path to the cache file for the isotope.

    """
    return _get_cache_dir() / f"{isotope_name}.csv"


def _parse_nuclide(isotope_name: str) -> str:
    """Parse an isotope name into the nuclide format expected by the IAEA database.

    Args:
        isotope_name: Name of the isotope to parse.
            Format should be 'ElementMass' (e.g. Ru106) or 'MassElement' (e.g. 106Ru).

    Returns:
        nuclide: Nuclide name in the format 'masselement' (e.g. '106ru').
        Note that the element symbol is converted to lowercase.

    """
    element, mass_number = _parse_isotope(isotope_name)
    return f"{mass_number}{element.lower()}"


def _download_spectrum_data(isotope_name: str) -> str:
    """Download the antineutrino spectrum for a given nuclide from the IAEA database.

    Args:
        isotope_name: Name of the isotope to download data for.
            Format should be 'ElementMass' (e.g. Ru106) or 'MassElement' (e.g. 106Ru).

    Returns:
        Path to the downloaded spectrum data file in the cache.

    """
    # Download the data file
    nuclide = _parse_nuclide(isotope_name)
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


def _load_spectrum_file(isotope_name: str) -> np.ndarray:
    """Load in antineutrino spectrum data from a CSV file in the cache.

    Args:
        isotope_name: Name of the isotope to load data for.
            Format should be 'ElementMass' (e.g. Ru106) or 'MassElement' (e.g. 106Ru).

    Returns:
        Array containing energy, flux, and uncertainty.

    """
    nuclide = _parse_nuclide(isotope_name)
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


def get_antineutrino_spectrum(isotope_name: str) -> np.ndarray:
    """Load in antineutrino spectrum data for a given isotope.

    If the spectrum data is not already in the cache, it is downloaded from the
    IAEA database and saved locally before loading.

    Args:
        isotope_name: Isotope name to load the spectrum for.
            Format should be 'ElementMass' (e.g. Ru106) or 'MassElement' (e.g. 106Ru).

    Returns:
        Array containing energy, flux, and uncertainty.

    """
    nuclide = _parse_nuclide(isotope_name)
    if not _get_cache_file(nuclide).is_file():
        _download_spectrum_data(nuclide)
    return _load_spectrum_file(nuclide)
