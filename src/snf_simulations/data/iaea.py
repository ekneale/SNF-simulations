"""Module for loading antineutrino spectrum data from the IAEA database."""

import os
import urllib.error
import urllib.request
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

from .utils import _parse_isotope

_CACHE_DIR_ENV_VAR = "SNF_SIMULATIONS_CACHE_DIR"


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
        "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0",
    )
    try:
        content = urllib.request.urlopen(req, timeout=5).read().decode("utf-8")  # noqa: S310
    except urllib.error.HTTPError as err:
        body = ""
        if err.fp is not None:
            body = err.fp.read().decode("utf-8", errors="ignore").lower()
        if err.code == 403 and "cloudflare" in body:  # noqa: PLR2004
            target_path = _get_cache_file(nuclide)
            msg = (
                f"IAEA request to url '{url}' was blocked by Cloudflare (HTTP 403).\n"
                "Try downloading the data from a web browser, and saving it to "
                f"the cache as '{target_path}'."
            )
            raise RuntimeError(msg) from err
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
