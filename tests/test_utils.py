"""Unit tests for utility functions."""

import numpy as np
import pytest

from snf_simulations.utils import linear_interpolate_with_errors, sample_histogram

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers


def test_linear_interpolate_with_errors_no_change() -> None:
    """Test interpolation returns original values on identical binning."""
    original_bins = np.array([0.0, 1.0, 2.0])
    original_content = np.array([10.0, 20.0])
    original_errors = np.array([1.0, 2.0])

    new_content, new_errors = linear_interpolate_with_errors(
        original_bins,
        original_content,
        original_errors,
        new_bins=original_bins,
    )

    assert np.allclose(new_content, original_content)
    assert np.allclose(new_errors, original_errors)


def test_linear_interpolate_with_errors_extrapolation() -> None:
    """Test bins fully outside original range are set to zero."""
    original_bins = np.array([0.0, 1.0, 2.0])
    original_content = np.array([10.0, 20.0])
    original_errors = np.array([1.0, 2.0])
    new_bins = np.array([-1.0, 0.0, 1.0, 2.0, 3.0])

    new_content, new_errors = linear_interpolate_with_errors(
        original_bins,
        original_content,
        original_errors,
        new_bins,
    )

    assert new_content[0] == 0.0
    assert new_content[-1] == 0.0
    assert new_errors[0] == 0.0
    assert new_errors[-1] == 0.0


def test_linear_interpolate_with_errors_overlap_uses_edge_error() -> None:
    """Test overlapping boundary bins use nearest edge-bin error."""
    original_bins = np.array([0.0, 1.0, 2.0])
    original_content = np.array([10.0, 20.0])
    original_errors = np.array([1.0, 2.0])
    new_bins = np.array([-0.2, 0.2, 1.2, 2.2])

    _, new_errors = linear_interpolate_with_errors(
        original_bins,
        original_content,
        original_errors,
        new_bins,
    )

    assert np.isclose(new_errors[0], original_errors[0])
    assert np.isclose(new_errors[-1], original_errors[-1])


def test_linear_interpolate_with_errors_invalid_inputs() -> None:
    """Test linear_interpolate_with_errors raises ValueError for invalid inputs."""
    with pytest.raises(
        ValueError,
        match=r"original_bins must have length len\(original_content\) \+ 1",
    ):
        _ = linear_interpolate_with_errors(
            original_bins=np.array([0.0, 1.0]),
            original_content=np.array([10.0, 20.0]),
            original_errors=np.array([1.0, 2.0]),
            new_bins=np.array([0.0, 1.0]),
        )

    with pytest.raises(
        ValueError,
        match="original_errors must have the same length as original_content",
    ):
        _ = linear_interpolate_with_errors(
            original_bins=np.array([0.0, 1.0, 2.0]),
            original_content=np.array([10.0, 20.0]),
            original_errors=np.array([1.0]),
            new_bins=np.array([0.0, 1.0]),
        )

    with pytest.raises(ValueError, match="original_bins must have at least two values"):
        _ = linear_interpolate_with_errors(
            original_bins=np.array([0.0]),
            original_content=np.array([]),  # need to be 1 less than original_bins
            original_errors=np.array([]),   # or else earlier check will fail first
            new_bins=np.array([0.0, 1.0]),
        )

    with pytest.raises(ValueError, match="new_bins must have at least two values"):
        _ = linear_interpolate_with_errors(
            original_bins=np.array([0.0, 1.0]),
            original_content=np.array([10.0]),
            original_errors=np.array([1.0]),
            new_bins=np.array([0.0]),
        )

    with pytest.raises(ValueError, match="original_bins must be strictly increasing"):
        _ = linear_interpolate_with_errors(
            original_bins=np.array([1.0, 0.0, 1.0]),
            original_content=np.array([10.0, 20.0]),
            original_errors=np.array([1.0, 2.0]),
            new_bins=np.array([0.0, 1.0]),
        )

    with pytest.raises(ValueError, match="new_bins must be strictly increasing"):
        _ = linear_interpolate_with_errors(
            original_bins=np.array([0.0, 1.0, 2.0]),
            original_content=np.array([10.0, 20.0]),
            original_errors=np.array([1.0, 2.0]),
            new_bins=np.array([1.0, 0.0, 1.0]),
        )


def test_sample_histogram_range() -> None:
    """Test sampling bounds compliance."""
    bin_edges = np.array([0.0, 1.0, 3.0])
    bin_contents = np.array([1.0, 2.0])

    samples = sample_histogram(bin_edges, bin_contents, samples=500)

    assert np.all(samples >= bin_edges[0])
    assert np.all(samples <= bin_edges[-1])


def test_sample_histogram_reproducible() -> None:
    """Test sampling reproducibility with seed."""
    bin_edges = np.array([0.0, 1.0, 3.0])
    bin_contents = np.array([1.0, 2.0])

    samples1 = sample_histogram(bin_edges, bin_contents, samples=500, seed=1234)
    samples2 = sample_histogram(bin_edges, bin_contents, samples=500, seed=1234)

    assert np.array_equal(samples1, samples2)


def test_sample_histogram_zero_area() -> None:
    """Test sampling fails for zero-weight histogram."""
    bin_edges = np.array([0.0, 1.0, 2.0])
    bin_contents = np.array([0.0, 0.0])

    with pytest.raises(
        ValueError, match="Histogram has zero total area; cannot sample"
    ):
        _ = sample_histogram(bin_edges, bin_contents, samples=10)


def test_sample_histogram_invalid_inputs() -> None:
    """Test sample_histogram raises ValueError for invalid inputs."""
    with pytest.raises(
        ValueError, match="bin_edges and bin_contents must be 1D arrays"
    ):
        _ = sample_histogram(
            np.array([[0.0, 1.0], [2.0, 3.0]]),
            np.array([1.0, 2.0]),
            samples=10,
        )

    with pytest.raises(
        ValueError, match="bin_edges and bin_contents must be 1D arrays"
    ):
        _ = sample_histogram(
            np.array([0.0, 1.0, 2.0]),
            np.array([[1.0, 2.0]]),
            samples=10,
        )

    with pytest.raises(
        ValueError,
        match=r"bin_edges must have length len\(bin_contents\) \+ 1",
    ):
        _ = sample_histogram(
            np.array([0.0, 1.0]),
            np.array([1.0, 2.0]),
            samples=10,
        )

    with pytest.raises(ValueError, match="bin_edges must be strictly increasing"):
        _ = sample_histogram(
            np.array([0.0, 0.0, 1.0]),
            np.array([1.0, 2.0]),
            samples=10,
        )

    with pytest.raises(ValueError, match="bin_contents must be non-negative"):
        _ = sample_histogram(
            np.array([0.0, 1.0, 2.0]),
            np.array([1.0, -1.0]),
            samples=10,
        )
