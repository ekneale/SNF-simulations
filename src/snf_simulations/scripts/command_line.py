"""Command line script to run SNF simulations and generate plots."""

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib.pyplot as plt

from snf_simulations.cask import Cask
from snf_simulations.data import get_example_tbq_path
from snf_simulations.detector import Detector
from snf_simulations.physics import calculate_flux_at_distance
from snf_simulations.spec import Spectrum


def print_detector_rates(spec: Spectrum, distance: float) -> None:
    """Calculate the antineutrino flux and expected event rates at a given distance."""
    # The total flux is the integral of the spectrum over the energy range of interest.
    # For inverse beta decay the threshold is 1.806 MeV, so we integrate above this
    # energy to get the total flux of antineutrinos that can be detected.
    total_flux = spec.integrate(lower_energy=1806)
    flux_at_distance = calculate_flux_at_distance(total_flux, distance=distance)
    print(f" Antineutrino flux at {distance} m:")
    print(f"  {flux_at_distance:.3e} per cm2 per second")
    print(f"  {flux_at_distance * 60 * 60 * 24:.3e} per cm2 per day")

    # The prototype detector is a plastic scintillator detector with a volume of
    # ~0.6 m^3, but we expect to use a minimum of 2 of these detectors to measure
    # a SNF cask so we'll start with 1.2 m^3 detector volume.
    volume = 0.6 * 2  # m^3
    proton_density = 4.6e22  # number density of protons in cm^-3
    detector = Detector(
        volume=volume,
        proton_density=proton_density,
        name="Baseline prototype",
    )

    # Get the expected event rate in the detector for this flux,
    # for lower and upper efficiencies.
    rate_lower = detector.calculate_event_rate(spec, distance, efficiency=0.3)
    rate_upper = detector.calculate_event_rate(spec, distance, efficiency=0.5)
    print(f" Event rate in {detector.name} detector:")
    print(f"  {rate_lower:.4e} to {rate_upper:.4e} per second")
    print(
        f"  {rate_lower * 60 * 60 * 24:.3f} to {rate_upper * 60**2 * 24:.3f} per day",
    )


def run_single(
    filepath: Path,
    cask_mass: float,
    simulation_times: Sequence[float],
    detector_distance: float = 40,
) -> dict[float, Spectrum]:
    """Plot the spectrum for a single dry cask of fuel after different cooling times.

    Args:
        filepath: Path to the input file.
        cask_mass: Total mass of the cask in kg.
        simulation_times: List of times to simulate the cask at, in years.
        detector_distance: Distance from the cask to the detector in meters.

    Returns:
        A dictionary mapping simulation times to Spectra for the single cask
            at each simulation time.

    """
    print("Generating single cask spectra after different cooling times...")

    # Create a single Cask of the given mass.
    cask = Cask.from_tabqfile(filepath, total_mass=cask_mass, name=filepath.stem)

    # Get the Spectra for the given times after removal from the core.
    spectra: dict[float, Spectrum] = {}
    for simulation_time in simulation_times:
        spec = cask.get_total_spectrum(cooling_time=simulation_time)
        spectra[simulation_time] = spec

    # Calculate and print flux and event rates for each simulation time.
    print(
        f"Single {cask_mass / 1000:.0f}-tonne {cask.name} cask at {detector_distance}m."
    )
    for simulation_time, spec in spectra.items():
        print(f"After {simulation_time} years:")
        print_detector_rates(spec, distance=detector_distance)

    # Save the first spectrum to CSV.
    print("Writing single cask spectrum data to CSV...")
    filename = f"{cask.name}_single.csv"
    spectra[simulation_times[0]].write_csv(output_filename=filename)
    print(f"Saved to {filename}")

    # Plot the cask spectra for each time.
    figure = plt.figure(figsize=(12, 6))
    axes = figure.add_subplot(1, 1, 1)
    for simulation_time, spec in spectra.items():
        spec.equalise(width=1, min_energy=0, max_energy=6000)
        axes.step(
            spec.energy[:-1],
            spec.flux,
            where="post",
            label=f"{simulation_time:.1f} years since removal from core",
        )
    axes.set_xlim(0, 6000)
    axes.set_ylim(1e10, 1e18)
    axes.set_xlabel("Energy [keV]")
    axes.set_ylabel("Relative Flux [keV^-1 s^-1]")
    axes.set_title(f"Single {cask_mass / 1000:.0f}-tonne {cask.name} cask spectrum")
    axes.set_yscale("log")
    axes.legend()

    print("Displaying plot...")
    plt.show()

    # Save the plot as a PDF.
    filename = f"{cask.name}_single.pdf"
    figure.savefig(filename)
    print(f"Saved to {filename}")

    return spectra


def run_multiple(
    filepath: Path,
    cask_mass: float,
    n_casks: int,
    cooling_times: Sequence[float],
    detector_distance: float = 40,
) -> Spectrum:
    """Calculate total spectra for multiple casks.

    Args:
        filepath: Path to the input file.
        cask_mass: Total mass of each cask in kg.
        n_casks: Number of casks to simulate.
        cooling_times: List of times the casks have been cooling, in years.
        detector_distance: Distance from the casks to the detector in meters.

    Returns:
        The total Spectrum for all the casks combined at the first cooling time.

    """
    print("Generating multiple cask spectra...")
    # Now we have multiple casks that have been cooling for different times,
    # and want to find the combined antineutrino flux from all of them.
    # Because each Cask has the same composition, and we have the same number per set,
    # we can multiply the total mass by the number of casks and just create
    # a single Cask object.
    cask = Cask.from_tabqfile(
        filepath, total_mass=cask_mass * n_casks, name=filepath.stem
    )

    # Create the Spectra for each set of casks at the specified times.
    spectra = []
    for cooling_time in cooling_times:
        spec = cask.get_total_spectrum(cooling_time=cooling_time)
        spectra.append(spec)

    # Combine the spectra from all casks to get the total spectrum for all 40 casks.
    spec_multiple = spectra[0]
    for spec in spectra[1:]:
        spec_multiple = spec_multiple + spec
    spec_multiple.name = cask.name

    # Save to CSV.
    print("Writing multiple cask spectrum data to CSV...")
    filename = f"{cask.name}_multiple.csv"
    spec_multiple.write_csv(output_filename=filename)
    print(f"Saved to {filename}")

    # Calculate and print flux and event rates for the combined spectrum.
    print(f"Multiple {cask.name} casks at {detector_distance}m:")
    print_detector_rates(spec_multiple, distance=detector_distance)

    return spec_multiple


def run_compare(
    spec_single: Spectrum,
    spec_multiple: Spectrum,
) -> None:
    """Plot both single and multiple flux spectra on one graph."""
    figure = plt.figure(figsize=(12, 6))
    axes = figure.add_subplot(1, 1, 1)

    spec_single.equalise(width=1, min_energy=0, max_energy=6000)
    axes.step(
        spec_single.energy[:-1],
        spec_single.flux,
        where="post",
        label=f"{spec_single.name}, single cask (0.5 yrs)",
    )

    spec_multiple.equalise(width=1, min_energy=0, max_energy=6000)
    axes.step(
        spec_multiple.energy[:-1],
        spec_multiple.flux,
        where="post",
        label=f"{spec_multiple.name}, multiple casks",
    )

    axes.set_xlim(0, 6000)
    axes.set_ylim(1e10, 1e18)
    axes.set_xlabel("Energy [keV]")
    axes.set_ylabel("Relative Flux [keV^-1 s^-1]")
    title = f"Comparison of single and multiple cask spectra for {spec_single.name}"
    axes.set_title(title)
    axes.set_yscale("log")
    axes.legend()

    print("Displaying plot...")
    plt.show()

    # Save the plot as a PDF.
    filename = f"{spec_single.name}_single.pdf"
    figure.savefig(filename)
    print(f"Saved to {filename}")


def run_multiple_full(
    filepath: Path,
    cask_mass: float,
    cask_cooling_times: Mapping[int, int] | Mapping[float, int],
    simulation_times: Sequence[float],
    detector_distance: float = 40,
) -> None:
    """Plot multiple flux spectra for different cooling times."""
    print("Generating multiple cask spectra for different cooling times...")
    # Now we do a full simulation of multiple sets of casks removed at different times,
    # and simulate their spectra into the future.
    # We don't have a set number of casks per set now, so we have to create
    # a Cask object for each set with the overall mass.
    casks = {}
    for cooling_time, n_casks in cask_cooling_times.items():
        casks[cooling_time] = Cask.from_tabqfile(
            filepath, total_mass=cask_mass * n_casks, name=filepath.stem
        )

    # Now get the spectra for all the sets of casks at each of the simulation times.
    spectra: dict[float, Spectrum] = {}
    for simulation_time in simulation_times:
        # We need to get the total spectrum for each set at this time,
        # and then combine them.
        time_spectra = []
        for cooling_time, cask in casks.items():
            new_cooling_time = cooling_time + simulation_time
            spec = cask.get_total_spectrum(cooling_time=new_cooling_time)
            time_spectra.append(spec)
        total_spec = time_spectra[0]
        for spec in time_spectra[1:]:
            total_spec = total_spec + spec
        total_spec.name = f"Total spectrum for all casks after {simulation_time} years"
        spectra[simulation_time] = total_spec

    # Calculate and print flux and event rates for each simulation time.
    print(
        f"Multiple {cask_mass / 1000:.0f}-tonne {cask.name} casks "
        f"at {detector_distance}m."
    )
    for simulation_time, spec in spectra.items():
        print(f"After {simulation_time} years:")
        print_detector_rates(spec, distance=detector_distance)

    # Plot the total spectra for each simulation time.
    figure = plt.figure(figsize=(12, 6))
    axes = figure.add_subplot(1, 1, 1)
    for simulation_time, spec in spectra.items():
        spec.equalise(width=1, min_energy=0, max_energy=6000)
        axes.step(
            spec.energy[:-1],
            spec.flux,
            where="post",
            label=f"Spectrum after {simulation_time} years",
        )
    axes.set_xlim(0, 6000)
    axes.set_ylim(1e10, 1e18)
    axes.set_xlabel("Energy [keV]")
    axes.set_ylabel("Relative Flux [keV^-1 s^-1]")
    axes.set_title(
        f"Multiple {cask_mass / 1000:.0f}-tonne {cask.name} cask spectra "
        f"at different simulation times",
    )
    axes.set_yscale("log")
    axes.legend()

    print("Displaying plot...")
    plt.show()

    # Save the plot as a PDF.
    filename = f"{cask.name}_multiple.pdf"
    figure.savefig(filename)
    print(f"Saved to {filename}")


def run_sample(spec: Spectrum, n_samples: int = 1000000) -> None:
    """Sample the cask spectrum to simulate detector observations."""
    print("Sampling cask spectrum...")

    # Take 1 million samples from the cask spectrum.
    spec.equalise(width=1, min_energy=0, max_energy=6000)
    samples = spec.sample(n_samples)

    # Plot the sampled spectra along with the original for comparison.
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
    axes.set_title(f"Comparison of original and sampled spectra for {spec.name}")
    axes.set_yscale("log")
    axes.legend()

    print("Displaying plot...")
    plt.show()

    # Save the plot as a PDF.
    filename = f"{spec.name}_sampled.pdf"
    figure.savefig(filename)
    print(f"Saved to {filename}")


def run(filepath: Path) -> None:
    """Run SNF simulations and calculate antineutrino fluxes and event rates.

    Performs calculations for single and multiple dry casks of spent nuclear fuel
    at various cooling times, calculating antineutrino flux and expected detector
    event rates at 40 meters distance.

    Args:
        filepath: Path to the input file.

    """
    # Run each of the test cases:
    # First, just a single cask at a set of different times after removal.
    cask_mass = 10000
    simulation_times = [0.5, 5, 10, 20]
    spec_single = run_single(filepath, cask_mass, simulation_times)

    # Then, multiple casks at the same removal times.
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    n_casks = 10
    spec_multiple = run_multiple(
        filepath,
        cask_mass,
        n_casks,
        simulation_times,
    )

    # Compare the two cases.
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    run_compare(spec_single[simulation_times[0]], spec_multiple)

    # Now do the full multiple cask case with sets of existing casks
    # hat have already been cooling for diferent times.
    # Here we have 10 casks that have just been removed (cooling_time=0),
    # 10 that have been cooling for 1 year, etc.
    # Then we simulate the total spectra for these casks
    # at each of the simulation times.
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    cask_cooling_times = {
        0: 10,
        1: 10,
        5: 10,
        10: 10,
        20: 10,
    }
    run_multiple_full(filepath, cask_mass, cask_cooling_times, simulation_times)

    # Finally, sample the spectrum to simulate detector observations.
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    run_sample(spec_multiple)


def main() -> None:
    """Parse command line arguments and run the simulation."""
    parser = argparse.ArgumentParser(
        description="Run SNF simulations and generate plots",
    )
    parser.add_argument(
        "filepath",
        type=Path,
        nargs="?",
        help="Path to the input file",
        default=get_example_tbq_path(),
    )
    args = parser.parse_args()
    run(filepath=args.filepath)


if __name__ == "__main__":
    main()
