# SNF-simulations
SNF_simulations simulates antineutrino emissions from dry casks of spent nuclear fuel by plotting a total energy spectrum dependent on the masses and cooling times of contributing isotopes commonly found in SNF.

## Installation
Use the package manager [pip](https://pip.pypa.io/en/stable/) to install SNF_simulations.

```bash
pip install SNF_simulations
```

## Usage
```python
import snf_simulations
```
Currently the package plots the antineutrino spectra of SNF originating from the Sizewell and Hartlepool reactors, both for one dry cask of SNF as well 15 dry casks of fuel with varying cooling times.

plotting.py combines most of the functions contained in the package. One output is the spectrum of one dry cask for either the Hartlepool or Sizewell reactor, depending on which is set to True, and how this spectrum varies over time. The other output is the total spectrum of the sum of 15 dry casks of SNF from either reactor. By running command_line.py these plots are created. The input proportions of isotopes can b=be varied in plotting.py, as well as removal times and number of dry casks contained at the storage facility. Isotopes can be added to the data directory, antineutrino_spec_data, and then also added into define_proportions.py. 
