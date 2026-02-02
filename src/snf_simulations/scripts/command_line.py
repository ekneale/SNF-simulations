"""Command line script to run SNF simulations and generate plots."""

from snf_simulations.physics import calculate_event_rate, calculate_flux
from snf_simulations.plotting import (
    multiple_fluxes,
    multiple_single_plot,
    plot_multiple,
    plot_multiple_casks_hartlepool,
    plot_multiple_casks_sizewell,
    plot_sample,
    plot_single_cask,
)
from snf_simulations.spec import write_spec


def main():
    """Run SNF simulations and generate plots."""
    # 1) Plot the spectrum for a single dry cask of fuel
    spec_single = plot_single_cask(
        reactor="sizewell",
        removal_times=[0.5, 1, 5, 10, 20],
    )
    data = write_spec(spec_single, output_filename="sizewell_single.csv")
    energy_single, flux_single = data[:, 0], data[:, 1]


    # 2) Plot multiple casks
    spec_multiple_sizewell = plot_multiple_casks_sizewell([0.5, 5, 10, 20])
    spec_multiple_hartlepool = plot_multiple_casks_hartlepool([3, 7, 15, 19])
    data = write_spec(
        spec_multiple_sizewell,
        output_filename="sizewell_multiple.csv",
    )  # change for sizewell or hartlepool
    energy_multiple, flux_multiple = data[:, 0], data[:, 1]


    # 3) Calculate flux at the detector located 40 m from the dry casks
    flux = calculate_flux(spec_multiple_sizewell, 40)
    print(flux, "cm^-2 s^-1")
    print(flux * 60 * 60 * 24, "cm^-2 days^-1")
    rate_lower, rate_upper = calculate_event_rate(flux, 0.2, 0.4)
    print(rate_lower, rate_upper, "per s")
    print(rate_lower * 60 * 60 * 24, rate_upper * 60**2 * 24, "per day")

    flux = calculate_flux(spec_single, 40)
    print(flux, "cm^-2 s^-1")
    print(flux * 60 * 60 * 24, "cm^-2 days^-1")
    rate_lower, rate_upper = calculate_event_rate(flux, 0.2, 0.4)
    print(rate_lower, rate_upper, "per s")
    print(rate_lower * 60 * 60 * 24, rate_upper * 60**2 * 24, "per day")


    # 4) Plot both flux spectra on one graph
    multiple_single_plot(
        energy_single,
        flux_single,
        energy_multiple,
        flux_multiple,
        reactor="sizewell",
    )

    # 5) Plot "lots_of_graphs" (TODO)
    lots_of_graphs = multiple_fluxes(reactor="sizewell")

    plot_multiple(lots_of_graphs, reactor="sizewell")

    for i in range(len(lots_of_graphs)):
        flux = calculate_flux(lots_of_graphs[i], 40)
        print(flux, "cm^-2 s^-1")
        print(flux * 60 * 60 * 24, "cm^-2 days^-1")
        rate_lower, rate_upper = calculate_event_rate(flux, 0.2, 0.4)
        print(rate_lower, rate_upper, "per s")
        print(rate_lower * 60 * 60 * 24, rate_upper * 60**2 * 24, "per day")

    # 6) Sample the multiple cask spectrum to simulate detector observations
    plot_sample(spec_multiple_sizewell, reactor="sizewell")


if __name__ == "__main__":
    main()
