"""Unit tests for data loading functions."""

from pathlib import Path

import numpy as np
import pytest

from snf_simulations.data import (
    _get_cache_dir,
    download_spectrum_data,
    get_reactors,
    load_antineutrino_data,
    load_isotope_data,
    load_reactor_data,
    load_spectrum,
    load_spectrum_file,
)

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers


def test_get_reactors() -> None:
    """Test that reactor list can be retrieved."""
    reactors = get_reactors()
    assert isinstance(reactors, list), "Reactors should be a list"
    assert len(reactors) > 0, "Should have at least one reactor"
    assert all(isinstance(r, str) for r in reactors), "All reactors should be strings"
    assert "sizewell" in reactors, "Sizewell should be in reactor list"
    assert "hartlepool" in reactors, "Hartlepool should be in reactor list"


def test_load_reactor_data_sizewell() -> None:
    """Test loading reactor data for Sizewell."""
    data = load_reactor_data("sizewell")
    assert isinstance(data, dict), "Loaded data should be a dictionary"
    assert len(data) > 0, "Should have loaded isotope data"
    assert all(isinstance(k, str) for k in data), "All keys should be isotope names"
    assert all(isinstance(v, float) for v in data.values()), (
        "All values should be floats"
    )
    assert all(v >= 0 for v in data.values()), "All proportions should be non-negative"


def test_load_reactor_data_hartlepool() -> None:
    """Test loading reactor data for Hartlepool."""
    data = load_reactor_data("hartlepool")
    assert isinstance(data, dict), "Loaded data should be a dictionary"
    assert len(data) > 0, "Should have loaded isotope data"


def test_load_reactor_data_invalid() -> None:
    """Test that loading invalid reactor raises ValueError."""
    with pytest.raises(ValueError, match=r"Reactor.*data file not found"):
        load_reactor_data("invalid_reactor")


def test_load_isotope_data_all() -> None:
    """Test loading all isotope data."""
    molar_masses, half_lives = load_isotope_data()
    assert isinstance(molar_masses, dict), "Molar masses should be a dictionary"
    assert isinstance(half_lives, dict), "Half lives should be a dictionary"
    assert len(molar_masses) > 0, "Should have loaded molar masses"
    assert len(half_lives) > 0, "Should have loaded half lives"
    assert len(molar_masses) == len(half_lives), (
        "Molar masses and half lives should have same length"
    )

    # Check data types and values
    assert all(isinstance(k, str) for k in molar_masses), (
        "All molar mass keys should be isotope names"
    )
    assert all(isinstance(v, int) for v in molar_masses.values()), (
        "All molar masses should be integers"
    )
    assert all(v > 0 for v in molar_masses.values()), (
        "All molar masses should be positive"
    )

    assert all(isinstance(k, str) for k in half_lives), (
        "All half life keys should be isotope names"
    )
    assert all(isinstance(v, float) for v in half_lives.values()), (
        "All half lives should be floats"
    )
    assert all(v > 0 for v in half_lives.values()), "All half lives should be positive"


def test_load_isotope_data() -> None:
    """Test loading specific isotopes."""
    isotopes = ["Sr90", "Y90"]
    molar_masses, half_lives = load_isotope_data(isotopes)
    assert len(molar_masses) == 2, "Should have loaded 2 isotopes"
    assert len(half_lives) == 2, "Should have loaded 2 isotopes"
    assert "Sr90" in molar_masses, "Sr90 should be in molar masses"
    assert "Y90" in molar_masses, "Y90 should be in molar masses"
    assert "Sr90" in half_lives, "Sr90 should be in half lives"
    assert "Y90" in half_lives, "Y90 should be in half lives"
    assert molar_masses["Sr90"] == 90, "Molar mass of Sr90 should be 90 g/mol"
    assert molar_masses["Y90"] == 90, "Molar mass of Y90 should be 90 g/mol"
    assert half_lives["Sr90"] == 28.91, (
        "Half life of Sr90 should be approximately 28.91 years"
    )
    assert half_lives["Y90"] == 0.0073, (
        "Half life of Y90 should be approximately 0.0073 years"
    )


def test_load_isotope_data_missing_file(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that missing isotope CSV raises ValueError with clear message."""

    # Point the data package root to an empty temporary directory so
    # isotopes.csv is absent.
    def _mock_files(_: str) -> Path:
        return tmp_path

    monkeypatch.setattr("snf_simulations.data.files", _mock_files)

    with pytest.raises(ValueError, match="Isotope CSV file not found."):
        load_isotope_data()


def test_load_spectrum() -> None:
    """Test loading spectrum data for a valid isotope."""
    spectrum = load_spectrum("Sr90")
    assert isinstance(spectrum, np.ndarray), "Spectrum should be a numpy array"
    assert spectrum.shape[1] == 3, (
        "Spectrum should have 3 columns (energy, flux, uncertainty)"
    )
    assert spectrum.shape[0] > 0, "Spectrum should have at least one row"

    # Check that energy values are positive and increasing
    energy = spectrum[:, 0]
    assert np.all(energy >= 0), "Energy values should be positive"
    assert np.all(np.diff(energy) > 0), "Energy should be increasing"

    # Check that flux values are non-negative
    flux = spectrum[:, 1]
    assert np.all(flux >= 0), "Flux values should be non-negative"

    # Check that uncertainties are non-negative
    uncertainty = spectrum[:, 2]
    assert np.all(uncertainty >= 0), "Uncertainty values should be non-negative"


def test_load_spectrum_invalid() -> None:
    """Test that loading invalid isotope raises ValueError."""
    with pytest.raises(ValueError, match=r"Nuclide format not recognized"):
        load_spectrum("InvalidIsotope")


def test_download_spectrum_data_uses_cache(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that downloaded spectra are written to a writable user cache."""
    csv_content = (
        "p_z,p_n,p_symbol,p_energy,d_z,d_n,d_symbol,bin_en,dn_de,unc_dn_de,"
        "dn_de_nu,unc_dn_de_nu,extraction_date\n"
        "38,52,Sr,0,39,51,Y,0.0,0.1,0.01,0.2,0.02,2026-04-15\n"
    )

    class _MockResponse:
        def read(self) -> bytes:
            return csv_content.encode("utf-8")

    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(
        "urllib.request.urlopen",
        lambda request, timeout=5: _MockResponse(),
    )

    filepath = Path(download_spectrum_data("Ru106"))

    assert filepath == tmp_path / "106ru.csv"
    assert filepath.is_file(), "Downloaded spectrum file should be written to cache"


def test_download_spectrum_data_wraps_network_errors(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that download errors are wrapped in RuntimeError."""
    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))

    def _raise_url_error(request: object, timeout: int = 5) -> object:
        del request, timeout
        raise TimeoutError("network down")

    monkeypatch.setattr("urllib.request.urlopen", _raise_url_error)

    with pytest.raises(RuntimeError, match=r"Error downloading spectrum data"):
        download_spectrum_data("Ru106")


def test_load_spectrum_file_missing_cache_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that loading a missing cached spectrum raises ValueError."""
    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))

    with pytest.raises(ValueError, match=r"not found in cache"):
        load_spectrum_file("Sr90")


def test_load_spectrum_downloads_when_cache_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that load_spectrum triggers a download when cache is missing."""
    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))
    called: list[str] = []

    def _fake_download(nuclide: str) -> str:
        called.append(nuclide)
        cache_file = tmp_path / f"{nuclide}.csv"
        cache_file.write_text(
            "p_z,p_n,p_symbol,p_energy,d_z,d_n,d_symbol,bin_en,dn_de,unc_dn_de,"
            "dn_de_nu,unc_dn_de_nu,extraction_date\n"
            "38,52,Sr,0,39,51,Y,0.5,0.1,0.01,1.5,0.15,2026-04-15\n",
            encoding="utf-8",
        )
        return str(cache_file)

    monkeypatch.setattr("snf_simulations.data.download_spectrum_data", _fake_download)

    spectrum = load_spectrum("Sr90")

    assert called == ["90sr"]
    np.testing.assert_allclose(spectrum, np.array([[0.5, 1.5, 0.15]]))


def test_get_cache_dir_uses_xdg_cache_home(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that the XDG cache location is used when no explicit cache is set."""
    monkeypatch.delenv("SNF_SIMULATIONS_CACHE_DIR", raising=False)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))

    cache_dir = _get_cache_dir()

    assert cache_dir == tmp_path / "snf_simulations" / "spec_data"
    assert cache_dir.is_dir()


def test_load_spectrum_uses_cached_file(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that cached spectra are loaded without attempting a download."""
    cache_file = tmp_path / "90sr.csv"
    cache_file.write_text(
        "p_z,p_n,p_symbol,p_energy,d_z,d_n,d_symbol,bin_en,dn_de,unc_dn_de,"
        "dn_de_nu,unc_dn_de_nu,extraction_date\n"
        "38,52,Sr,0,39,51,Y,1.5,0.1,0.01,2.5,0.25,2026-04-15\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(
        "snf_simulations.data.download_spectrum_data",
        lambda nuclide: pytest.fail(f"Unexpected download for {nuclide}"),
    )

    spectrum = load_spectrum("Sr90")

    np.testing.assert_allclose(spectrum, np.array([[1.5, 2.5, 0.25]]))


def test_load_antineutrino_data() -> None:
    """Test loading antineutrino data for multiple isotopes."""
    isotopes = ["Sr90", "Y90"]
    data = load_antineutrino_data(isotopes)
    assert isinstance(data, dict), "Loaded data should be a dictionary"
    assert len(data) == 2, "Should have loaded 2 isotopes"

    for isotope in isotopes:
        assert isotope in data, f"{isotope} should be in loaded data"
        spectrum = data[isotope]
        assert isinstance(spectrum, np.ndarray), (
            f"{isotope} spectrum should be numpy array"
        )
        assert spectrum.shape[1] == 3, f"{isotope} spectrum should have 3 columns"
