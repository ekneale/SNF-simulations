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

Some notes on using this code:
As previously mentioned command_line.py runs the rest of the code. If you want to get out the energy/flux csv files for a single cask from either reactor and whatever removal time, line 5 of command_line.py should not be hashed out, if you want the energy/flux csv files for multiple casks from either reactor and whatever respective removal time, line 12 should not be hashed out. Make sure you only have one of these lines running at one time.
It's also important that the reactor you have set to true remains consisitent throughout command_line.py so make sure of that.
In flux.py, within the write_spec function which writes the energy/flux csvs, change the file name on line 47 to either single cask or multiple cask and also the number for whatever removal time you are looking at.
The reactor name should select automatically so you don't need to change that manually.
In plotting.py, within the plot_single_cask function, make sure to return the spectra.At() with the right number for the removal time you are focussing on. For example if you were looking at 10 years, you would return spectra.At(3) as it's the third item in the list. To see the lists of removal times, refer to command_line.py. [0.5,1,5,10,20]
This is similar in the plot_multiple_casks_sizewell and plot_multiple_casks_hartlepool functions. You need to change the return casks.At() for whatever removal time you want. The removal time lists are different for sizewell [0.5,5,10,20] and hartlepool [3,7,15,19] for the multiple casks so refer to command_line.py where you should see these lists I've given and make sure to select the correct item (removal time). 
