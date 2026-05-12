# SNF-simulations
SNF-simulations generates antineutrino spectra that can be used to simulate the antineutrino emission from spent nuclear fuel (SNF).

For more information about what this package does and how to use it, see the [documentation](https://ekneale.github.io/SNF-simulations/index.html).

## Installation
To use, clone the repository then use the package manager [pip](https://pip.pypa.io/en/stable/) to install SNF_simulations.

```bash
pip install snf-simulations
```

If you want to run the dashboard locally, you will need to install the extra dependencies using the following command:

```bash
pip install snf-simulations[dashboard]
```

### Dependencies

SNF-simulations depends on the following packages:

- `numpy`
- `pandas`
- `mendeleev` (for accessing periodic table data)
- `matplotlib` (for plotting with the `snf-sim` demo script)

All dependencies are automatically installed when you install SNF-simulations with pip.

With the `dashboard` option, the following packages are also installed:

- `shiny` (the dashboard is built using the [Shiny framework for Python](https://shiny.posit.co/py/))
- `shinywidgets` (for interactive widgets in the dashboard)
- `plotly` (for interactive plots)
