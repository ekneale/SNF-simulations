"""Data loading module for antineutrino spectra calculations."""

from .fispin import get_isotope_masses
from .iaea import get_antineutrino_spectrum
from .mendeleev import get_isotope_properties
from .utils import get_example_tbq_path

__all__ = [
    "get_antineutrino_spectrum",
    "get_example_tbq_path",
    "get_isotope_masses",
    "get_isotope_properties",
]
