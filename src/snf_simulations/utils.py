"""Utility functions for spectrum interpolation and sampling."""

import numpy as np


def linear_interpolate_with_errors(
    original_bins: np.ndarray,
    original_content: np.ndarray,
    original_errors: np.ndarray,
    new_bins: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Linearly interpolate histogram content and propagate errors onto new bins."""
    if len(original_bins) != len(original_content) + 1:
        msg = "original_bins must have length len(original_content) + 1"
        raise ValueError(msg)
    if len(original_errors) != len(original_content):
        msg = "original_errors must have the same length as original_content"
        raise ValueError(msg)
    if len(original_bins) < 2:  # noqa: PLR2004
        msg = "original_bins must have at least two values"
        raise ValueError(msg)
    if len(new_bins) < 2:  # noqa: PLR2004
        msg = "new_bins must have at least two values"
        raise ValueError(msg)
    if np.any(np.diff(original_bins) <= 0):
        msg = "original_bins must be strictly increasing"
        raise ValueError(msg)
    if np.any(np.diff(new_bins) <= 0):
        msg = "new_bins must be strictly increasing"
        raise ValueError(msg)

    # Interpolate new content values
    original_centres = (original_bins[:-1] + original_bins[1:]) / 2
    new_centres = (new_bins[:-1] + new_bins[1:]) / 2
    new_content = np.interp(new_centres, original_centres, original_content)

    # Any extrapolated bins outside the original range should be set to zero
    lower_edges = new_bins[:-1]
    upper_edges = new_bins[1:]
    lower_mask = upper_edges <= original_bins[0]
    upper_mask = lower_edges >= original_bins[-1]
    extrapolation_mask = lower_mask | upper_mask
    new_content[extrapolation_mask] = 0

    # Propagate errors
    new_errors = np.zeros_like(new_centres)
    for i, centre in enumerate(new_centres):
        # Check for extrapolated bins - these should keep zero error
        if extrapolation_mask[i]:
            continue

        # Handle boundary centres explicitly to avoid out-of-range indexing.
        # For bins that overlap the original edge range but whose centres fall outside
        # the original centre range, keep the nearest edge-bin error constant.
        # This avoids pseudo-bin (under/overflow) interpolation at the boundaries.
        if centre <= original_centres[0] or np.isclose(centre, original_centres[0]):
            new_errors[i] = original_errors[0]
            continue
        if centre >= original_centres[-1] or np.isclose(centre, original_centres[-1]):
            new_errors[i] = original_errors[-1]
            continue

        # Find the two closest original bin centers.
        # Using np.searchsorted finds the "insertion point" for the new centre,
        # i.e. the index of where it would go to keep the array sorted.
        # So if the original_centres are [1, 2, 3] and centre is 2.5, idx will be 2
        # as it would fit between 2 (index 1) and 3 (index 2).
        # Therefore the surrounding bins are at idx-1 and idx.
        idx = int(np.searchsorted(original_centres, centre))
        lower_idx = idx - 1
        upper_idx = idx

        # Calculate new error by propagating errors from the two surrounding bins,
        # weighted by distance to the new centre.
        c_lower = original_centres[lower_idx]
        c_upper = original_centres[upper_idx]
        err_lower = original_errors[lower_idx]
        err_upper = original_errors[upper_idx]
        weight_lower = (c_upper - centre) / (c_upper - c_lower)
        weight_upper = (centre - c_lower) / (c_upper - c_lower)
        new_errors[i] = np.sqrt(
            weight_lower**2 * err_lower**2 + weight_upper**2 * err_upper**2,
        )

    return new_content, new_errors


def sample_histogram(
    bin_edges: np.ndarray,
    bin_contents: np.ndarray,
    samples: int = 100,
    seed: int | None = None,
) -> np.ndarray:
    """Sample x values from histogram bins, similar to ROOT TH1::GetRandom.

    Args:
        bin_edges: 1D array of bin edges with length N+1.
        bin_contents: 1D array of bin contents with length N.
        samples: Number of samples to draw.
        seed: Seed for reproducible random sampling.

    Returns:
        Array of sampled x values.

    """
    if bin_edges.ndim != 1 or bin_contents.ndim != 1:
        msg = "bin_edges and bin_contents must be 1D arrays"
        raise ValueError(msg)
    if len(bin_edges) != len(bin_contents) + 1:
        msg = "bin_edges must have length len(bin_contents) + 1"
        raise ValueError(msg)
    widths = np.diff(bin_edges)
    if np.any(widths <= 0):
        msg = "bin_edges must be strictly increasing"
        raise ValueError(msg)
    if np.any(bin_contents < 0):
        msg = "bin_contents must be non-negative"
        raise ValueError(msg)

    # Match ROOT TH1::GetRandom behaviour: bin selection probability is
    # proportional to bin content, then sample uniformly within the selected bin.
    weights = bin_contents
    total_weight = np.sum(weights)
    if total_weight <= 0:  # Avoid division by zero errors
        msg = "Histogram has zero total area; cannot sample"
        raise ValueError(msg)
    probabilities = weights / total_weight

    # Use numpy's random choice to select X bins according to their probabilities,
    # for the requested number of samples.
    rng = np.random.default_rng(seed)
    sampled_indices = rng.choice(len(bin_contents), size=samples, p=probabilities)

    # Finally, for each bin take a uniform sample between the upper and lower edges.
    # This gives a continuous distribution of sampled x values from within the bins.
    lower = bin_edges[sampled_indices]
    upper = bin_edges[sampled_indices + 1]
    return rng.uniform(lower, upper)
