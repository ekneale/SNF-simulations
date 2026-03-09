"""Command line script to run SNF simulations and generate plots."""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from snf_simulations.cask import Cask
from snf_simulations.data import load_reactor_data
from snf_simulations.physics import calculate_event_rate, calculate_flux_at_distance
from snf_simulations.spec import Spectrum


def run_single(reactor: str = "sizewell", cask_mass: float = 10000) -> Spectrum:
    """Plot the spectrum for a single dry cask of fuel at different removal times."""
    print("Generating single cask spectra at different cooling times...")

    # Create a single Cask of the given mass,
    # using the isotope proportions for the given reactor
    isotope_proportions = load_reactor_data(reactor)
    cask = Cask(
        isotope_proportions=isotope_proportions,
        total_mass=cask_mass,
        name=f"{reactor}_cask",
    )

    # Create the Spectra for removal times of 0, 0.5, 1, 5, 10, 20 years.
    removal_times = [0, 0.5, 1, 5, 10, 20]
    spectra = []
    for removal_time in removal_times:
        spec = cask.get_total_spectrum(removal_time=removal_time)
        spectra.append(spec)

    # Calculate and print flux and event rates for the 0.5 year removal time
    spec_05 = spectra[1]
    # The total flux is the integral of the spectrum over the energy range of interest.
    # For inverse beta decay the threshold is 1.806 MeV, so we integrate above this
    # energy to get the total flux of antineutrinos that can be detected.
    total_flux = spec_05.integrate(lower_energy=1806)
    flux_at_40m = calculate_flux_at_distance(total_flux, distance=40)
    print()
    print(
        f"Single {cask_mass / 1000:.0f}-tonne {reactor.capitalize()} cask "
        "at 40m after 0.5 years:"
    )
    print(" Antineutrino flux:")
    print(f"  {flux_at_40m:.3e} per cm2 per second")
    print(f"  {flux_at_40m * 60 * 60 * 24:.3e} per cm2 per day")
    rate_lower, rate_upper = calculate_event_rate(flux_at_40m, 0.2, 0.4)
    print(" Event rate in VIDARR detector:")
    print(f"  {rate_lower:.3e} to {rate_upper:.3e} per second")
    print(
        f"  {rate_lower * 60 * 60 * 24:.3f} to {rate_upper * 60**2 * 24:.3f} per day",
    )

    # Save 0.5 year spectrum to CSV
    print()
    print("Writing single cask spectrum data to CSV...")
    filename = f"{reactor.capitalize()}_single.csv"
    spec.write_csv(output_filename=filename)
    print(f"Saved to {filename}")

    # Plot the cask spectra after different removal times
    figure = plt.figure(figsize=(12, 6))
    axes = figure.add_subplot(1, 1, 1)
    for spec, removal_time in zip(spectra, removal_times, strict=True):
        spec.equalise(width=1, min_energy=0, max_energy=6000)
        axes.step(
            spec.energy[:-1],
            spec.flux,
            where="post",
            label=f"{removal_time:.1f} years since removal from core",
        )
    axes.set_xlim(0, 6000)
    axes.set_ylim(1e10, 1e18)
    axes.set_xlabel("Energy [keV]")
    axes.set_ylabel("Relative Flux [keV^-1 s^-1]")
    axes.set_title(
        f"Single {cask_mass / 1000:.0f}-tonne {reactor.capitalize()} cask spectrum"
    )
    axes.set_yscale("log")
    axes.legend()

    print()
    print("Displaying plot...")
    plt.show()

    # Save the plot as a PDF
    filename = f"{reactor.capitalize()}_single.pdf"
    figure.savefig(filename)
    print(f"Saved to {filename}")

    return spec_05


def run_multiple(
    reactor: str = "sizewell",
    cask_mass: float = 10000,
    casks_per_removal: int = 10,
) -> Spectrum:
    """Calculate total spectra for multiple casks."""
    print("Generating multiple cask spectra...")
    # 40 casks in total, each 10 tonnes, split between 4 different removal times.
    # Since we're combining the spectra at the end, instead of ten 10-tonne casks
    # at each removal time, we just simulate one 100-tonne cask.
    # Since all the casks are the same mass and proportions, just different removal
    # times, we only need to create a single Cask object.
    isotope_proportions = load_reactor_data(reactor)
    cask = Cask(
        isotope_proportions=isotope_proportions,
        total_mass=cask_mass * casks_per_removal,
        name=f"{reactor}_cask",
    )

    # Create the Spectra for each cask at each removal time
    if reactor == "sizewell":
        removal_times = [0.5, 5, 10, 20]
    elif reactor == "hartlepool":
        removal_times = [3, 7, 15, 19]
    spectra = []
    for removal_time in removal_times:
        spec = cask.get_total_spectrum(removal_time=removal_time)
        spectra.append(spec)

    # Combine the spectra from all casks to get the total spectrum for all 40 casks.
    spec_multiple = spectra[0]
    for spec in spectra[1:]:
        spec_multiple = spec_multiple + spec
    spec_multiple.name = "Total spectrum for all casks"

    # Save to CSV
    print()
    print("Writing multiple cask spectrum data to CSV...")
    filename = f"{reactor.capitalize()}_multiple.csv"
    spec_multiple.write_csv(output_filename=filename)
    print(f"Saved to {filename}")

    # Calculate and print flux and event rates for the combined spectrum
    total_flux = spec_multiple.integrate(lower_energy=1806)
    flux_at_40m = calculate_flux_at_distance(total_flux, distance=40)
    print()
    print(f"Multiple {reactor.capitalize()} casks at 40m:")
    print(" Antineutrino flux:")
    print(f"  {flux_at_40m:.3e} per cm2 per second")
    print(f"  {flux_at_40m * 60 * 60 * 24:.3e} per cm2 per day")
    rate_lower, rate_upper = calculate_event_rate(flux_at_40m, 0.2, 0.4)
    print(" Event rate in VIDARR detector:")
    print(f"  {rate_lower:.3e} to {rate_upper:.3e} per second")
    print(
        f"  {rate_lower * 60 * 60 * 24:.3f} to {rate_upper * 60**2 * 24:.3f} per day",
    )

    return spec_multiple


def run_compare(
    spec_05: Spectrum,
    spec_multiple: Spectrum,
    reactor: str = "sizewell",
) -> None:
    """Plot both single and multiple flux spectra on one graph."""
    figure = plt.figure(figsize=(12, 6))
    axes = figure.add_subplot(1, 1, 1)

    spec_05.equalise(width=1, min_energy=0, max_energy=6000)
    axes.step(
        spec_05.energy[:-1],
        spec_05.flux,
        where="post",
        label=f"{reactor.capitalize()}, single cask (0.5 yrs)",
    )

    spec_multiple.equalise(width=1, min_energy=0, max_energy=6000)
    axes.step(
        spec_multiple.energy[:-1],
        spec_multiple.flux,
        where="post",
        label=f"{reactor.capitalize()}, multiple casks",
    )

    axes.set_xlim(0, 6000)
    axes.set_ylim(1e10, 1e18)
    axes.set_xlabel("Energy [keV]")
    axes.set_ylabel("Relative Flux [keV^-1 s^-1]")
    axes.set_title(
        f"Comparison of single and multiple cask spectra for {reactor.capitalize()}"
    )
    axes.set_yscale("log")
    axes.legend()

    print()
    print("Displaying plot...")
    plt.show()

    # Save the plot as a PDF
    filename = f"{reactor.capitalize()}_single.pdf"
    figure.savefig(filename)
    print(f"Saved to {filename}")


def run_multiple_full(
    reactor: str = "sizewell",
    cask_mass: float = 10000,
    casks_per_removal: int = 10,
) -> None:
    """Plot multiple flux spectra for different cooling times."""
    print("Generating multiple cask spectra for different cooling times...")
    # This is the same setup as run_multi, but now we add extra cooling times
    # to see how the spectrum evolves.
    isotope_proportions = load_reactor_data(reactor)
    cask = Cask(
        isotope_proportions=isotope_proportions,
        total_mass=cask_mass * casks_per_removal,
        name=f"{reactor}_cask",
    )

    if reactor == "sizewell":
        removal_times = [0.5, 5, 10, 20]
    elif reactor == "hartlepool":
        removal_times = [3, 7, 15, 19]
    cooling_times = [0, 1, 5, 10, 20]

    # Create the combined spectra from all casks for each cooling time
    spectra = []
    for cooling_time in cooling_times:
        new_ages = np.array(removal_times) + cooling_time
        time_spectra = []
        for removal_time in new_ages:
            spec = cask.get_total_spectrum(removal_time=removal_time)
            time_spectra.append(spec)
        total_spec = time_spectra[0]
        for spec in time_spectra[1:]:
            total_spec = total_spec + spec
        total_spec.name = f"Total spectrum for all casks after {cooling_time} years"
        spectra.append(total_spec)

    # Calculate and print flux and event rates for each cooling time
    print()
    for spec, cooling_time in zip(spectra, cooling_times, strict=True):
        total_flux = spec.integrate(lower_energy=1806)
        flux_at_40m = calculate_flux_at_distance(total_flux, distance=40)
        print(
            f"Multiple cask flux at 40 m for {reactor.capitalize()}"
            f" after {cooling_time} years:",
        )
        print(" Antineutrino flux:")
        print(f"  {flux_at_40m:.3e} per cm2 per second")
        print(f"  {flux_at_40m * 60 * 60 * 24:.3e} per cm2 per day")
        rate_lower, rate_upper = calculate_event_rate(flux_at_40m, 0.2, 0.4)
        print(" Event rate in VIDARR detector:")
        print(f"  {rate_lower:.3e} to {rate_upper:.3e} per second")
        print(
            f"  {rate_lower * 60 * 60 * 24:.3f} to "
            f"{rate_upper * 60**2 * 24:.3f} per day",
        )
        print()

    # Plot the total spectra for each cooling time
    figure = plt.figure(figsize=(12, 6))
    axes = figure.add_subplot(1, 1, 1)
    for spec, cooling_time in zip(spectra, cooling_times, strict=True):
        spec.equalise(width=1, min_energy=0, max_energy=6000)
        axes.step(
            spec.energy[:-1],
            spec.flux,
            where="post",
            label=f"Spectrum after {cooling_time} years",
        )
    axes.set_xlim(0, 6000)
    axes.set_ylim(1e10, 1e18)
    axes.set_xlabel("Energy [keV]")
    axes.set_ylabel("Relative Flux [keV^-1 s^-1]")
    axes.set_title(
        f"Multiple {cask_mass / 1000:.0f}-tonne {reactor.capitalize()} cask spectra "
        f"at different cooling times",
    )
    axes.set_yscale("log")
    axes.legend()

    print()
    print("Displaying plot...")
    plt.show()

    # Save the plot as a PDF
    filename = f"{reactor.capitalize()}_multiple.pdf"
    figure.savefig(filename)
    print(f"Saved to {filename}")


def run_sample(
    spec: Spectrum,
    reactor: str = "sizewell",
) -> None:
    """Sample the cask spectrum to simulate detector observations."""
    print("Sampling cask spectrum...")

    # Take 1 million samples from the cask spectrum
    spec.equalise(width=1, min_energy=0, max_energy=6000)
    samples = spec.sample(samples=1000000)

    # Plot the sampled spectra along with the original for comparison
    figure = plt.figure(figsize=(12, 6))
    axes = figure.add_subplot(1, 1, 1)

    samples_hist = axes.hist(
        samples,
        bins=list(spec.energy[:-1]),
        label=f"Sampled spectrum ({len(samples)} samples)",
        histtype="step",
        align="mid",
    )

    scaled_spec = spec * (max(samples_hist[0]) / max(spec.flux))
    axes.step(
        scaled_spec.energy[:-1],
        scaled_spec.flux,
        where="post",
        label="Original spectrum (scaled)",
        ls="dashed",
    )

    axes.set_xlim(0, 6000)
    axes.set_xlabel("Energy [keV]")
    axes.set_ylabel("Relative Flux [keV^-1 s^-1]")
    axes.set_title(
        f"Comparison of original and sampled spectra for {reactor.capitalize()}"
    )
    axes.set_yscale("log")
    axes.legend()

    print()
    print("Displaying plot...")
    plt.show()

    # Save the plot as a PDF
    filename = f"{reactor.capitalize()}_sampled.pdf"
    figure.savefig(filename)
    print(f"Saved to {filename}")


def run(reactor: str = "sizewell", cask_mass: float = 10000) -> None:
    """Run SNF simulations and calculate antineutrino fluxes and event rates.

    Performs calculations for single and multiple dry casks of spent nuclear fuel
    at various cooling times, calculating antineutrino flux and expected detector
    event rates at 40 meters distance.

    Args:
        reactor: Reactor name ('sizewell' or 'hartlepool').
        cask_mass: Total mass of SNF in each cask (kg).

    """
    # Run each of the test cases
    spec_single = run_single(reactor=reactor, cask_mass=cask_mass)
    print()
    spec_multiple = run_multiple(reactor=reactor, cask_mass=cask_mass)
    print()
    run_compare(spec_single, spec_multiple, reactor=reactor)
    print()
    run_multiple_full(reactor=reactor, cask_mass=cask_mass)
    print()
    run_sample(spec_multiple, reactor=reactor)


def main() -> None:
    """Parse command line arguments and run the simulation."""
    parser = argparse.ArgumentParser(
        description="Run SNF simulations and generate plots",
    )
    parser.add_argument(
        "--reactor",
        type=str,
        choices=["sizewell", "hartlepool"],
        default="sizewell",
        help="Reactor name: 'sizewell' or 'hartlepool' (default: sizewell)",
    )
    args = parser.parse_args()
    run(reactor=args.reactor.lower())


if __name__ == "__main__":
    main()
