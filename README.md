# SNF-simulations
This package simulates antineutrino emissions from dry casks of spent nuclear fuel (SNF). This is done by plotting a total energy spectrum dependent on the masses and cooling times of isotopes commonly found in SNF and a calculation of the expected flux of antineutrinos at the detector.

## Installation
To use, clone the repository then use the package manager [pip](https://pip.pypa.io/en/stable/) to install SNF_simulations.

```bash
pip install -i https://test.pypi.org/simple/ snf-simulations
```

## Usage
```python
import snf_simulations
```
Currently, the antineutrino_spec_data directory contains data for the antineutrino emission spectrum for 16 isotopes found in SNF. By running command_line.py several spectra are plotted. The first is the antineutrino spectrum varying with cooling times for one dry cask of SNF originating from either the Hartlepool detector or Sizewell detecor. The second is the total spectrum for 10 dry casks x 10 tonnes x 4 cooling times, so the simulated spectrum at an interim storage facility containing 40 dry casks of SNF. The plot following this shows that spectrum varying with time assuming no dry casks are added. The final spectrum is the sampled total spectrum for the simulated observations of an interim storage facility. A simple and approximate calculation of flux and event rate is also performed.

The modules can be used to create antineutrino emission spectra for a single dry cask as well as multiple. The module plotting.py provides examples on how to do this, by combining all modules. Spectra for different reactors can be plotted by chaning input mass proportions for the isotopes contained in the antineutrino_spec_data diredtory. 

More isotopes can be added to the antineutrino_spec_data directory. In this case the define_proportions module would have to be updated, as well as load_data to incldue additional isotopes. 
