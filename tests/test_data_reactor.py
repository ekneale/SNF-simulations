"""Unit tests for reactor data functions."""

import pytest

from snf_simulations.data.reactor import get_reactor_data

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers


def test_get_reactor_data_sizewell() -> None:
    """Test loading reactor data for Sizewell."""
    data = get_reactor_data("sizewell")
    assert isinstance(data, dict), "Loaded data should be a dictionary"
    assert len(data) > 0, "Should have loaded isotope data"
    assert all(isinstance(k, str) for k in data), "All keys should be isotope names"
    assert all(isinstance(v, float) for v in data.values()), (
        "All values should be floats"
    )
    assert all(v >= 0 for v in data.values()), "All proportions should be non-negative"


def test_get_reactor_data_hartlepool() -> None:
    """Test loading reactor data for Hartlepool."""
    data = get_reactor_data("hartlepool")
    assert isinstance(data, dict), "Loaded data should be a dictionary"
    assert len(data) > 0, "Should have loaded isotope data"


def test_get_reactor_data_invalid() -> None:
    """Test that loading invalid reactor raises ValueError."""
    with pytest.raises(ValueError, match=r"Reactor.*data file not found"):
        get_reactor_data("invalid_reactor")
