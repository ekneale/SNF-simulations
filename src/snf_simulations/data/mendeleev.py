"""Module for loading isotope data from the mendeleev package."""

from mendeleev import isotope

from .utils import _UNITS_TO_SECONDS, _parse_isotope


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
    half_life_years = half_life * seconds_per_unit / _UNITS_TO_SECONDS["year"]

    return {"molar_mass": molar_mass, "half_life": half_life_years}
