# SNF-simulations
SNF_simulations simulates antineutrino emissions from dry casks of spent nuclear fuel by plotting a total energy spectrum dependent on the mass proportions of isotopes found in the SNF and cooling time.

## Installation
Use the package manager [pip](https://pip.pypa.io/en/stable/) to install SNF_simulations.

```bash
pip install SNF_simulations
```

## Usage
```python
import snf_simulations

# returns antineutrino energy spectrum data for 16 isotpoes contained in antineutrino_spec_data directory
snf_simulations.load_data.load_data()

# returns antineutrino energy spectrum plotted as a histogram for a selected isotope
snf_simulations.load_spec.load_spec(E,dN,errors, "isotope")

#returns antineutrino energy spectrum plotted as a histogram with equal bin widths of 1 keV for a selected isotope
snf_simulations.load_spec.load_equal("name","isotope", E, dN, error, max_E, min_E=-0)
