"""Unit tests for sampling spectra."""

import os

import numpy as np
import ROOT

from snf_simulations.sample import sample_spec

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts


def _make_spec():
    spec = ROOT.TH1D("spec", "", 3, 0, 3)
    spec.SetBinContent(1, 1.0)
    spec.SetBinContent(2, 2.0)
    spec.SetBinContent(3, 3.0)
    return spec


def test_sample_spec_returns_samples():
    """Test sample size and bounds from a ROOT spectrum."""
    ROOT.gRandom.SetSeed(1234)  # Set seed for reproducibility
    spec = _make_spec()
    samples = sample_spec(spec, samples=5)

    # Samples should be fixed given the seed
    samples_ref = [
        1.0745583504904062,
        1.9929909987840801,
        2.244217532686889,
        2.6356768854893744,
        1.81318321172148,
    ]
    assert samples == samples_ref, "Sampled values do not match reference"


def test_sample_spec_writes_csv():
    """Test CSV output writing when a filename is provided."""
    ROOT.gRandom.SetSeed(1234)  # Set seed for reproducibility
    spec = _make_spec()
    filename = "samples"
    samples = sample_spec(spec, samples=3, output_filename=str(filename))

    # Check that the file was created and contents are correct
    # (note that .csv extension is added automatically by sample_spec)
    csv_file = filename + ".csv"
    saved = np.loadtxt(str(csv_file), delimiter=",")
    assert saved.tolist() == samples, "Saved samples do not match returned samples"

    # Clean up the file after testing
    os.remove(csv_file)
