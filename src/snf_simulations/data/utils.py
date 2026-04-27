"""Utility functions for loading and processing data."""

import re


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
