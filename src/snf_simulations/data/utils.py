"""Utility functions for loading and processing data."""

import re

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
