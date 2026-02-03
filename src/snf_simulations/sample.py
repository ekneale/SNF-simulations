"""Module for sampling a spectrum histogram to simulate detector observations."""

import numpy as np
import ROOT


def sample_spec(
    spec: ROOT.TH1D,
    samples: int = 100,
    output_filename: str | None = None,
) -> list[float]:
    """Sample the total spectrum to simulate what a detector could observe.

    Samples can be saved as CSV for further analysis.

    Args:
        spec: The flux spectrum histogram.
        samples: Number of samples to draw from the spectrum.
        output_filename: Filename to save the sampled data as a CSV file.
            If None, data is not saved.

    Returns:
        List of sampled flux values.

    """
    sampled_flux = [spec.GetRandom() for _ in range(samples)]

    if output_filename:
        # Save samples as a CSV file
        if not output_filename.endswith(".csv"):
            output_filename += ".csv"
        np.savetxt(
            output_filename,
            sampled_flux,
            delimiter=",",
            comments="",
        )

    return sampled_flux
