---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
# The `Spectrum` Class

## Introduction

The {py:obj}`Spectrum <snf_simulations.spec.Spectrum>` class, defined in the {py:obj}`snf_simulations.spec <snf_simulations.spec>` module, is the core data class used in the package to represent the antineutrino spectrum emitted by spent nuclear fuel (SNF).

Each `Spectrum` is a histogram, with energy bins as the x-axis and the corresponding antineutrino flux values, and uncertainties, as the y-axis. This class provides methods for manipulating and analyzing the spectrum data, such as equalising, integration and sampling.


## Creating a `Spectrum`

You can create a new `Spectrum` object by providing the necessary parameters: energy bins (in keV), flux values (in keV^-1) and uncertainties on the flux values. Each of these should be Numpy arrays.

:::{important}
When giving an array of energy values, the values are defined as the edges of the histogram bins, meaning that if you have $N$ bins, you should provide $N+1$ energy values.
:::

For example:

```{code-cell}python3
import numpy as np
from snf_simulations.spec import Spectrum

# Input arrays
energy_bins = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
flux_values = np.array([10, 20, 15, 5])
flux_errors = np.array([1, 2, 1.5, 0.5])

# Create a Spectrum object
spectrum = Spectrum(energy_bins, flux_values, flux_errors)
print(spectrum)
```

### Creating a `Spectrum` for a known isotope

In practice, often be easier to create a `Spectrum` object for a known isotope using the {py:func}`from_isotope <snf_simulations.spec.Spectrum.from_isotope>` class method, which retrieves the spectrum data from the internal database. For example:

```{code-cell}python3
spectrum = Spectrum.from_isotope("Am242")
print(spectrum)
```

This uses the {py:func}`load_spectrum <snf_simulations.data.load_spectrum>` function to load the spectrum data from one of the included antineutrino spectrum data files. For isotopes that are not included in the database, this will raise a `ValueError`, and you will need to create the `Spectrum` object manually using the constructor as shown in the previous example.

:::{note}
In the future, isotope data will be downloaded directly from the IAEA. See [this GitHub Issue](https://github.com/ekneale/SNF-simulations/issues/16) for details.
:::


## Manipulating `Spectrum` Objects

The `Spectrum` class provides several methods for manipulating the spectrum data.

### Equalising a `Spectrum`

As Spectrum objects are histograms, they can be equalised to a different set of energy bins using the {py:func}`equalise <snf_simulations.spec.Spectrum.equalise>` method. This is useful when you want to compare spectra that have different binning, or when you want to re-bin a spectrum to match the binning of another spectrum.

The {py:func}`equalise <snf_simulations.spec.Spectrum.equalise>` method takes three optional parameters: `width`, `min_energy` and `max_energy`. If `width` is provided, the method will create new energy bins with the specified width, starting from `min_energy` and ending at `max_energy`. If `width` is not provided, the method will default to bins of width 1 keV. `min_energy` and `max_energy` default to the minimum and maximum energy values of the original spectrum, respectively.

For example:

```{code-cell}python3
energy_bins = np.array([0, 2, 4, 6, 8])
flux_values = np.array([10, 20, 15, 5])
flux_errors = np.array([1, 2, 1.5, 0.5])

spectrum = Spectrum(energy_bins, flux_values, flux_errors)

# Equalise the spectrum to new bins
spectrum.equalise(width=1, max_energy=10)
print(spectrum.energy)
print(spectrum.flux)
print(spectrum.errors)
```

Note that any bins outside of the original energy range will be filled with zeros, and the corresponding uncertainties will also be set to zero.

### Scaling a `Spectrum`

Spectrum objects can be scaled by a constant factor using the `*` operator. This can be useful for normalising spectra, or to scale by an activity level. For example:

```{code-cell}python3
energy_bins = np.array([0, 2, 4, 6, 8])
flux_values = np.array([10, 20, 15, 5])
flux_errors = np.array([1, 2, 1.5, 0.5])

spectrum = Spectrum(energy_bins, flux_values, flux_errors)

# Scale the spectrum by a factor of 2
scaled_spectrum = spectrum * 2
print(scaled_spectrum.flux)
print(scaled_spectrum.errors)
```


### Integrating a `Spectrum`

The {py:func}`integrate <snf_simulations.spec.Spectrum.integrate>` method can be used to calculate the total flux in a specified energy range. This is done by integrating the flux values over the specified energy range, taking into account the bin widths and uncertainties. The method takes two parameters: `lower_energy` and `upper_energy`, which define the energy range for the integration. For example:

```{code-cell}python3
energy_bins = np.array([0, 2, 4, 6, 8])
flux_values = np.array([10, 20, 15, 5])
flux_errors = np.array([1, 2, 1.5, 0.5])

spectrum = Spectrum(energy_bins, flux_values, flux_errors)

# Integrate the spectrum from 1 keV to 5 keV
total_flux = spectrum.integrate(lower_energy=1, upper_energy=5)
print(f"Total flux: {total_flux}")
```


### Adding `Spectrum` objects

Two `Spectrum` objects can be added together using the `+` operator, which will add the flux values and combine the uncertainties in quadrature. For example:

```{code-cell}python3
energy_bins1 = np.array([0, 2, 4, 6, 8])
flux_values1 = np.array([10, 20, 15, 5])
flux_errors1 = np.array([1, 2, 1.5, 0.5])

spectrum1 = Spectrum(energy_bins1, flux_values1, flux_errors1)

energy_bins2 = np.array([0, 2, 4, 6, 8])
flux_values2 = np.array([5, 10, 15, 20])
flux_errors2 = np.array([0.5, 1, 1.5, 2])

spectrum2 = Spectrum(energy_bins2, flux_values2, flux_errors2)

# Add the two spectra together
combined_spectrum = spectrum1 + spectrum2
print(combined_spectrum.flux)
print(combined_spectrum.errors)
```

The energy bins of the two spectra **must** be the same for the addition to work. If they are not, you can use the {py:func}`equalise <snf_simulations.spec.Spectrum.equalise>` method to re-bin one of the spectra to match the other before adding them together.


### Sampling from a `Spectrum`

The {py:func}`sample <snf_simulations.spec.Spectrum.sample>` method can be used to generate random samples of antineutrino energies based on the flux values in the spectrum. This is done by treating the flux values as a probability distribution and sampling from it accordingly. The method takes one parameter: `samples`, which specifies the number of random samples to generate. For example:

```{code-cell}python3
energy_bins = np.array([0, 2, 4, 6, 8])
flux_values = np.array([10, 20, 15, 5])
flux_errors = np.array([1, 2, 1.5, 0.5])

spectrum = Spectrum(energy_bins, flux_values, flux_errors)

# Generate 10 random samples
samples = spectrum.sample(samples=10)
print(samples)
```

### Saving and Loading `Spectrum` objects

`Spectrum` objects can be saved to a file using the {py:func}`save <snf_simulations.spec.Spectrum.write_csv>` method, which saves the spectrum data in a simple text format. The method takes one parameter: `filename`, which specifies the name of the file to save the spectrum to. For example:

```{code-cell}python3
energy_bins = np.array([0, 2, 4, 6, 8])
flux_values = np.array([10, 20, 15, 5])
flux_errors = np.array([1, 2, 1.5, 0.5])

spectrum = Spectrum(energy_bins, flux_values, flux_errors)

# Save the spectrum to a file
spectrum.write_csv("spectrum.csv")
```

`Spectrum` objects can be loaded from a file using the {py:func}`load <snf_simulations.spec.Spectrum.from_file>` class method, which reads the spectrum data from a file in the same format as the `save` method. The method takes one parameter: `filename`, which specifies the name of the file to load the spectrum from. For example:

```{code-cell}python3
# Load a spectrum from a file
loaded_spectrum = Spectrum.from_file("spectrum.csv")
print(loaded_spectrum)
```
