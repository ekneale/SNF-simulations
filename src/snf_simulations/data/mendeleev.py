"""Module for loading isotope data from the mendeleev package."""

from functools import cache
from typing import TypedDict

import numpy as np
from mendeleev import isotope

from .utils import _UNITS_TO_SECONDS, _parse_isotope


class IsotopeProperties(TypedDict):
    """Class to represent the properties of an isotope."""

    molar_mass: float
    half_life: float
    decay_modes: list[str]


@cache
def _get_isotope_properties_cached(isotope_name: str) -> IsotopeProperties:
    """Return cached isotope properties loaded from mendeleev."""
    element, mass_number = _parse_isotope(isotope_name)
    mendeleev_isotope = isotope(element, mass_number)

    molar_mass = float(mendeleev_isotope.mass)  # ty: ignore
    if mendeleev_isotope.half_life is not None:
        half_life = float(mendeleev_isotope.half_life)  # ty: ignore
        unit = str(mendeleev_isotope.half_life_unit)
        seconds_per_unit = _UNITS_TO_SECONDS.get(unit)
        if seconds_per_unit is None:
            msg = f"Unsupported half-life unit for isotope {isotope_name}: {unit}. "
            msg += "Supported units are: "
            msg += ", ".join(sorted(_UNITS_TO_SECONDS))
            raise ValueError(msg)
        half_life_years = half_life * seconds_per_unit / _UNITS_TO_SECONDS["year"]
    else:
        half_life_years = np.inf
    decay_modes = [d.mode for d in mendeleev_isotope.decay_modes]

    return IsotopeProperties(
        molar_mass=molar_mass, half_life=half_life_years, decay_modes=decay_modes
    )


def get_isotope_properties(isotope_name: str) -> IsotopeProperties:
    """Get the mass, half-life and decay modes for the given isotope.

    Uses data from the mendeleev package.

    Args:
        isotope_name: Name of the isotope.
            Format should be 'ElementMass' (e.g. Ru106) or 'MassElement' (e.g. 106Ru).

    Returns:
        Dictionary containing:
        - the molar mass of the isotope (in g/mol)
        - the half-life of the isotope (in years)
        - the decay modes of the isotope (as a list of strings)

    """
    return _get_isotope_properties_cached(isotope_name)
