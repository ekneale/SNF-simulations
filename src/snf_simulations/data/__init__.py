"""Data loading module for antineutrino spectra calculations."""

from .fispin import get_isotope_proportions
from .iaea import get_antineutrino_spectrum
from .mendeleev import get_isotope_properties

__all__ = [
    "get_antineutrino_spectrum",
    "get_isotope_properties",
    "get_isotope_proportions",
]
