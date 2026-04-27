"""Module for loading isotope data from the mendeleev package."""

from mendeleev import isotope

from .utils import _parse_isotope

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
