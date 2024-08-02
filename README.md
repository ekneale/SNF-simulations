# SNF-simulations
This package simulates antineutrino emissiond from dry casks of spent nuclear fuel (SNF). This is done by plotting a total energy spectrum dependent on the masses and cooling times of isotopes commonly found in SNF and a calculation of the expected flux of antineutrinos at the detector.

## Installation
To use, clone the repository then use the package manager [pip](https://pip.pypa.io/en/stable/) to install SNF_simulations.

```bash
pip install SNF_simulations
```

## Usage
```python
import snf_simulations
```
Currently, the antineutrino_spec_data directory contains data for the antineutrino emission spectrum for 16 isotopes found in SNF. By running command_line.py two spectra are plotted The first is the antineutrino spectrum varying with cooling times for one dry cask of SNF originating from the Hartlepool detector, however, within the module if Sizewell is set to true the resulting spectrum will represent antineutrino emissions from a dry cask of fuel originiatiing from teh SIzewell reactor. The second spectrum displays the total antineutrino spectrum produced by 15 Hartlepool dry casks of SNF with different times since removal from the core, similarly this can also be plotted for Sizewell. Additionally, a calculation of expected flux at the detector is performed.

The modules can be used to create antineutrino emission spectra for dry casks originating from different reactors. This is done using  define_proportions.TotSpec(), where the inputs include the removal time for the dry cask of SNF and proportions of isotopes contained in the database, which are by default set to zero. Multiple casks can be added using 



plotting.py combines most of the functions contained in the package. One output is the spectrum of one dry cask for either the Hartlepool or Sizewell reactor, depending on which is set to True, and how this spectrum varies over time. The other output is the total spectrum of the sum of 15 dry casks of SNF from either reactor. By running command_line.py these plots are created. The input proportions of isotopes can b=be varied in plotting.py, as well as removal times and number of dry casks contained at the storage facility. Isotopes can be added to the data directory, antineutrino_spec_data, and then also added into define_proportions.py. 
