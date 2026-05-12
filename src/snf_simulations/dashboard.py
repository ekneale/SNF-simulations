"""Interactive Shiny dashboard for SNF antineutrino spectrum visualization."""

from pathlib import Path
from typing import TypedDict

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shinywidgets import output_widget, render_widget

from snf_simulations.cask import Cask
from snf_simulations.data import get_example_tbq_path
from snf_simulations.physics import calculate_event_rate, calculate_flux_at_distance


class SimInputs(TypedDict):
    """Snapshot of sidebar inputs for recalculated cask simulations."""

    filepath: Path
    cask_name: str
    cask_mass: float
    n_casks: int
    cooling_times: list[float]
    detector_distance: float


# Here we define the UI layout of the dashboard using Shiny's components.
# The sidebar holds the simulation inputs, and the main page has tabs for different
# output plots.
app_ui = ui.page_fluid(
    ui.busy_indicators.use(spinners=True, pulse=True, fade=True),
    ui.busy_indicators.options(
        spinner_delay="150ms",
        spinner_type="ring",
        spinner_size="2rem",
    ),
    ui.head_content(
        ui.tags.title("SNF Antineutrino Spectrum Dashboard"),
        ui.tags.style(
            # Custom CSS to increase spacing in the sidebar and reduce clipping,
            # and to style the header and footer.
            """
            .bslib-sidebar-layout > button.collapse-toggle {
                display: none !important;
            }
            .bslib-sidebar-layout > .sidebar > .sidebar-content {
                padding-top: 1rem !important;
            }
            .sidebar-content.bslib-gap-spacing {
                gap: 16px !important;
            }
            .sidebar-content.bslib-gap-spacing > h4 {
                margin-bottom: 0 !important;
            }
            .sidebar-content.bslib-gap-spacing >
                .form-group.shiny-input-container:first-of-type {
                margin-top: -4px;
                margin-bottom: -4px;
            }
            #tbq_file_progress {
                height: auto;
                overflow: visible;
                margin-top: 0.5rem;
            }
            #cooling_times.shiny-input-container > label {
                margin-bottom: 0.75rem;
            }
            .snf-header {
                display: flex;
                align-items: center;
                justify-content: space-between;
                padding: 0.75rem 1.25rem;
                background-color: #2c3e50;
                color: white;
            }
            .snf-header h1 {
                font-size: 1.4rem;
                margin: 0;
                font-weight: 600;
            }
            .snf-header-links {
                display: flex;
                gap: 1rem;
                align-items: center;
            }
            .snf-header-links a {
                color: rgba(255,255,255,0.85);
                text-decoration: none;
                font-size: 0.9rem;
            }
            .snf-header-links a:hover {
                color: white;
                text-decoration: underline;
            }
            .snf-footer {
                padding: 0.6rem 1.25rem;
                background-color: #2c3e50;
                color: rgba(255,255,255,0.75);
                font-size: 0.8rem;
                flex-shrink: 0;
            }

            """
        ),
    ),
    ui.tags.header(
        {"class": "snf-header"},
        ui.tags.h1("SNF Antineutrino Spectrum Dashboard"),
        ui.tags.div(
            {"class": "snf-header-links"},
            ui.tags.a(
                "📚 Documentation",
                href="https://ekneale.github.io/SNF-simulations",
                target="_blank",
                rel="noopener noreferrer",
            ),
            ui.tags.a(
                "📦 GitHub",
                href="https://github.com/ekneale/SNF-simulations",
                target="_blank",
                rel="noopener noreferrer",
            ),
        ),
    ),
    ui.layout_sidebar(
        ui.sidebar(
            ui.h4("Input Settings"),
            ui.input_file(
                "tbq_file",
                "Upload .tbQ file:",
                accept=[".tbQ"],
            ),
            ui.input_numeric(
                "cask_mass",
                "Cask mass (kg):",
                value=10000,
                min=1,
                step=1000,
            ),
            ui.input_slider(
                "n_casks",
                "Number of casks:",
                min=1,
                max=50,
                value=1,
            ),
            ui.input_checkbox_group(
                "cooling_times",
                "Cooling times (years):",
                {"0.5": "0.5", "1": "1", "5": "5", "10": "10"},
                selected=["0.5", "1", "5", "10"],
            ),
            ui.input_text(
                "new_cooling_time",
                "Add cooling time (years):",
                placeholder="e.g. 2.5",
            ),
            ui.input_action_button("add_cooling_time", "Add Cooling Time"),
            ui.input_slider(
                "detector_distance",
                "Detector distance (m):",
                min=1,
                max=200,
                value=40,
                step=1,
            ),
            ui.input_action_button("recalculate", "Recalculate", class_="btn-success"),
            width="250px",
        ),
        ui.navset_tab(
            ui.nav_panel(
                "Cask Simulations",
                output_widget("plot_cask_simulations"),
                ui.output_table("table_cask_rates"),
            ),
            ui.nav_panel(
                "Component Spectra",
                ui.input_select(
                    "component_cooling_time",
                    "Cooling time (years):",
                    choices={"0.5": "0.5", "1": "1", "5": "5", "10": "10"},
                    selected="0.5",
                ),
                output_widget("plot_component_spectra"),
            ),
            ui.nav_panel(
                "Spectrum Sampling",
                ui.input_select(
                    "sampling_cooling_time",
                    "Cooling time (years):",
                    choices={"0.5": "0.5", "1": "1", "5": "5", "10": "10"},
                    selected="0.5",
                ),
                ui.input_slider(
                    "n_samples",
                    "Number of samples:",
                    min=10000,
                    max=1000000,
                    value=100000,
                    step=10000,
                ),
                output_widget("plot_sampling"),
            ),
        ),
    ),
    ui.tags.footer(
        {"class": "snf-footer"},
        "\u00a9 Copyright 2026, E Kneale.",
    ),
)


def server(input: Inputs, output: Outputs, session: Session) -> None:  # noqa: PLR0915
    """Server logic for the Shiny dashboard."""
    # -----------------------------------------------------------------------
    # Simulation inputs
    sim_inputs: reactive.Value[SimInputs] = reactive.value(
        # Reactive value to hold the current simulation parameters
        {
            "filepath": get_example_tbq_path(),
            "cask_name": get_example_tbq_path().stem,
            "cask_mass": 10000.0,
            "n_casks": 1,
            "cooling_times": [0.5, 1.0, 5.0, 10.0],
            "detector_distance": 40.0,
        }
    )

    @reactive.calc
    def cask_reactive() -> Cask | None:
        """Reactive function to create a Cask instance from the current sim inputs."""
        try:
            params = sim_inputs()
            filepath = params["filepath"]
            cask = Cask.from_tabqfile(
                filepath,
                total_mass=float(params["cask_mass"]),
                name=params["cask_name"],
            )
            return cask
        except Exception as e:
            ui.notification_show(
                f"Error loading .tbQ file: {e}",
                type="error",
                duration=5,
            )
            return None

    @reactive.effect
    @reactive.event(input.recalculate)
    def _update_sim_inputs() -> None:
        """Capture sidebar values when the recalculate button is clicked."""
        try:
            if input.tbq_file() is None:
                filepath = get_example_tbq_path()
                cask_name = filepath.stem
            else:
                file_info = input.tbq_file()
                filepath = Path(file_info[0]["datapath"])
                cask_name = Path(file_info[0]["name"]).stem
            cask_mass = float(input.cask_mass())
            n_casks = int(input.n_casks())
            cooling_times = [float(t) for t in input.cooling_times()]
            detector_distance = float(input.detector_distance())

            snapshot: SimInputs = {
                "filepath": filepath,
                "cask_name": cask_name,
                "cask_mass": cask_mass,
                "n_casks": n_casks,
                "cooling_times": cooling_times,
                "detector_distance": detector_distance,
            }
            sim_inputs.set(snapshot)
        except Exception as e:
            ui.notification_show(
                f"Error recalculating cask simulations: {e}",
                type="error",
                duration=5,
            )

    # -----------------------------------------------------------------------
    # Cooling time options
    cooling_time_options: reactive.Value[list[float]] = reactive.value(
        # Reactive value to hold the list of cooling time options in the side bar
        [0.5, 1.0, 5.0, 10.0]
    )

    @reactive.effect
    @reactive.event(input.add_cooling_time)
    def _add_cooling_time() -> None:
        """Add a new cooling time option from the sidebar text input."""
        # Check the input isn't blank
        raw_value = (input.new_cooling_time() or "").strip()
        if not raw_value:
            ui.notification_show(
                "Please enter a cooling time value before adding.",
                type="warning",
                duration=3,
            )
            return

        # Check the input is a valid float
        try:
            new_value = float(raw_value)
        except ValueError:
            ui.notification_show(
                f"Invalid cooling time: {raw_value}",
                type="error",
                duration=4,
            )
            return

        # Check the input is non-negative
        if new_value < 0:
            ui.notification_show(
                "Cooling time must be non-negative.",
                type="error",
                duration=4,
            )
            return

        # Check the input is not earlier than the current cask age
        cask = cask_reactive()
        if cask is not None and new_value < cask.initial_cooling_time:
            ui.notification_show(
                (
                    "Cooling time was below the current cask age, so it was "
                    f"set to {cask.initial_cooling_time:.3e} years "
                    f"({cask.initial_cooling_time * 365.25:.1f} days)."
                ),
                type="warning",
                duration=4,
            )
            new_value = cask.initial_cooling_time

        # Check the input isn't already in the list
        current_options = cooling_time_options()
        tolerance = 1e-3  # Tolerance for duplicates
        if any(abs(value - new_value) < tolerance for value in current_options):
            ui.notification_show(
                f"Cooling time {new_value:g} is already in the list.",
                type="message",
                duration=3,
            )
            ui.update_text("new_cooling_time", value="")
            return

        # Add the new value and update the checkbox options (note we sort first)
        updated_options = sorted([*current_options, new_value])
        cooling_time_options.set(updated_options)

        # Now update the displayed choices, keeping the selected boxes the same
        # (note the new option is always selected)
        options = [f"{value:g}" for value in updated_options]
        selected_options = list(input.cooling_times() or [])
        if f"{new_value:g}" not in selected_options:
            selected_options.append(f"{new_value:g}")
        ui.update_checkbox_group(
            "cooling_times",
            choices={value: value for value in options},
            selected=[value for value in selected_options if value in options],
        )

        # Clear the text input
        ui.update_text("new_cooling_time", value="")

        # Update component cooling time dropdown
        selected = (
            input.component_cooling_time()
            if input.component_cooling_time() in options
            else options[0]
        )
        ui.update_select(
            "component_cooling_time",
            choices={value: value for value in options},
            selected=selected,
        )

        # Update sampling cooling time dropdown
        selected = (
            input.sampling_cooling_time()
            if input.sampling_cooling_time() in options
            else options[0]
        )
        ui.update_select(
            "sampling_cooling_time",
            choices={value: value for value in options},
            selected=selected,
        )

    # -----------------------------------------------------------------------
    # Cask Simulations tab
    @output
    @render_widget
    def plot_cask_simulations() -> go.Figure | None:
        """Plot cask spectra at selected cooling times for selected cask count."""
        cask = cask_reactive()
        if cask is None:
            return None

        params = sim_inputs()
        cask_mass = float(params["cask_mass"])
        n_casks = int(params["n_casks"])
        cooling_times = list(params["cooling_times"])

        # Create Plotly figure
        figure = go.Figure()
        figure.update_layout(
            title=f"{n_casks} × {cask_mass / 1000:.0f}-tonne {cask.name} cask spectra",
            xaxis={
                "title": "Energy [keV]",
                "range": [0, 6000],
                "unifiedhovertitle": {"text": "Energy: %{x:.0f} keV"},
            },
            yaxis={
                "title": "Relative Flux [keV⁻¹ s⁻¹]",
                "type": "log",
                "range": [10, 18],
            },
            template="plotly_white",
            hovermode="x unified",
        )

        for cooling_time in cooling_times:
            # Get the spectrum for this cooling time,
            # scale by the number of casks and equalise
            spec = cask.get_total_spectrum(cooling_time=cooling_time)
            spec = spec * n_casks
            spec.equalise(width=1, min_energy=0, max_energy=6000)

            # Add to the plot
            figure.add_trace(
                go.Scatter(
                    x=spec.energy[:-1],
                    y=spec.flux,
                    mode="lines",
                    line={"shape": "hv"},
                    name=f"{cooling_time:.1f} years since removal",
                    hovertemplate=(
                        f"Cooling time: {cooling_time:.1f} years<br>"
                        "Flux: %{y:.3e} keV⁻¹ s⁻¹"
                        "<extra></extra>"
                    ),
                )
            )

        figure.update_layout(legend={"title": {"text": "Cooling time"}})

        return figure

    @output
    @render.table
    def table_cask_rates() -> pd.DataFrame:
        """Table of flux and event rates for selected cask count and cooling times."""
        cask = cask_reactive()
        if cask is None:
            return pd.DataFrame()

        params = sim_inputs()
        n_casks = int(params["n_casks"])
        cooling_times = list(params["cooling_times"])
        detector_distance = float(params["detector_distance"])

        rows = []

        for cooling_time in cooling_times:
            spec = cask.get_total_spectrum(cooling_time=float(cooling_time))
            spec = spec * n_casks
            total_flux = spec.integrate(lower_energy=1806)
            flux_at_distance = calculate_flux_at_distance(
                total_flux, distance=detector_distance
            )
            rate_lower, rate_upper = calculate_event_rate(flux_at_distance, 0.2, 0.4)

            rows.append(
                {
                    "Cooling Time (yrs)": cooling_time,
                    "Flux (cm⁻²s⁻¹)": f"{flux_at_distance:.2e}",
                    "Flux (cm⁻²day⁻¹)": f"{flux_at_distance * 86400:.2e}",
                    "Event Rate Lower (s⁻¹)": f"{rate_lower:.2e}",
                    "Event Rate Upper (s⁻¹)": f"{rate_upper:.2e}",
                    "Event Rate Lower (day⁻¹)": f"{rate_lower * 86400:.2e}",
                    "Event Rate Upper (day⁻¹)": f"{rate_upper * 86400:.2e}",
                }
            )

        return pd.DataFrame(rows)

    # -----------------------------------------------------------------------
    # Component Spectra tab
    @output
    @render_widget
    def plot_component_spectra() -> go.Figure | None:
        """Plot per-isotope component spectra."""
        cask = cask_reactive()
        if cask is None:
            return None

        # Inputs
        cooling_time = float(input.component_cooling_time() or "0.5")

        # Get component spectra for this cooling time
        component_specs = cask.get_component_spectra(cooling_time=cooling_time)

        # Create Plotly figure
        figure = go.Figure()
        figure.update_layout(
            xaxis={
                "title": "Energy [keV]",
                "range": [0, 6000],
                "unifiedhovertitle": {"text": "Energy: %{x:.0f} keV"},
            },
            yaxis={
                "title": "Relative Flux [keV⁻¹ s⁻¹]",
                "type": "log",
                "range": [8, 18],
            },
            template="plotly_white",
            hovermode="closest",
            legend={"title": {"text": "Isotope"}},
        )

        # Plot each isotope spectrum
        for spec in component_specs:
            spec.equalise(width=1, min_energy=0, max_energy=6000)
            figure.add_trace(
                go.Scatter(
                    x=spec.energy[:-1],
                    y=spec.flux,
                    mode="lines",
                    line={"shape": "hv"},
                    name=spec.name,
                    opacity=0.7,
                    hovertemplate=(
                        "Isotope: %{fullData.name}<br>"
                        "Energy: %{x:.0f} keV<br>"
                        "Flux: %{y:.3e} keV⁻¹ s⁻¹"
                        "<extra></extra>"
                    ),
                )
            )
        return figure

    # -----------------------------------------------------------------------
    # Spectrum Sampling tab
    @output
    @render_widget
    def plot_sampling() -> go.Figure | None:
        """Plot sampled spectrum vs original."""
        cask = cask_reactive()
        if cask is None:
            return None

        # Inputs
        cooling_time = float(input.sampling_cooling_time() or "0.5")
        n_samples = input.n_samples()

        # Get spectrum and sample
        spec = cask.get_total_spectrum(cooling_time=cooling_time)
        spec.equalise(width=1, min_energy=0, max_energy=6000)
        samples = spec.sample(n_samples=n_samples)
        counts, _ = np.histogram(samples, bins=spec.energy)
        bin_centres = 0.5 * (spec.energy[:-1] + spec.energy[1:])

        # Create Plotly figure
        figure = go.Figure()
        figure.update_layout(
            xaxis={
                "title": "Energy [keV]",
                "range": [0, 6000],
                "unifiedhovertitle": {"text": "Energy: %{x:.0f} keV"},
            },
            yaxis={
                "title": "Counts / scaled flux",
                "type": "log",
            },
            template="plotly_white",
            hovermode="x unified",
        )

        # Plot sampled spectrum
        figure.add_trace(
            go.Scatter(
                x=bin_centres,
                y=np.maximum(counts, 1),
                mode="lines",
                line={"shape": "hv"},
                name="Sampled spectrum",
                hovertemplate=("Counts: %{y:.0f}<extra></extra>"),
            )
        )

        # Plot original spectrum for comparison
        scale_factor = max(counts.max(), 1) / max(spec.flux.max(), 1)
        scaled_flux = np.maximum(spec.flux * scale_factor, 1)
        figure.add_trace(
            go.Scatter(
                x=spec.energy[:-1],
                y=scaled_flux,
                mode="lines",
                line={"shape": "hv", "dash": "dash", "width": 2},
                name="Original (scaled)",
                hovertemplate=("Scaled flux: %{y:.3e}<extra></extra>"),
            )
        )

        return figure


app = App(app_ui, server)
