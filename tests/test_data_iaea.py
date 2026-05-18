"""Unit tests for IAEA antineutrino spectrum data functions."""

import urllib.error
from io import BytesIO
from pathlib import Path

import numpy as np
import pytest

from snf_simulations.data.iaea import (
    _copy_packaged_spectrum_to_cache,
    _download_spectrum_data,
    _get_cache_dir,
    _load_spectrum_file,
    get_antineutrino_spectrum,
)

# Suppress assert warnings from ruff
# ruff: noqa: S101  # asserts
# ruff: noqa: PLR2004  # magic numbers


def test_get_cache_dir_uses_xdg_cache_home(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that the XDG cache location is used when no explicit cache is set."""
    monkeypatch.delenv("SNF_SIMULATIONS_CACHE_DIR", raising=False)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))

    cache_dir = _get_cache_dir()

    assert cache_dir == tmp_path / "snf_simulations" / "spec_data"
    assert cache_dir.is_dir()


def test_get_cache_dir_defaults_to_home(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test cache defaults to ~/.cache when env vars are not set."""
    monkeypatch.delenv("SNF_SIMULATIONS_CACHE_DIR", raising=False)
    monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
    monkeypatch.setattr("snf_simulations.data.iaea.Path.home", lambda: tmp_path)

    cache_dir = _get_cache_dir()

    assert cache_dir == tmp_path / ".cache" / "snf_simulations" / "spec_data"
    assert cache_dir.is_dir()


def test_copy_packaged_spectrum_to_cache_copies_existing_file(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test bundled spectrum files are copied into the writable cache."""
    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))

    copied = _copy_packaged_spectrum_to_cache("Sr90")

    assert copied
    copied_file = tmp_path / "90sr.csv"
    assert copied_file.is_file()
    header = copied_file.read_text(encoding="utf-8").splitlines()[0]
    assert "bin_en" in header
    assert "dn_de_nu" in header


def test_copy_packaged_spectrum_to_cache_returns_false_when_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test helper returns False when there is no bundled file for isotope."""
    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))

    copied = _copy_packaged_spectrum_to_cache("Xe999")

    assert not copied
    assert not (tmp_path / "999xe.csv").exists()


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

    filepath = Path(_download_spectrum_data("Ru106"))

    assert filepath == tmp_path / "106ru.csv"
    assert filepath.is_file(), "Downloaded spectrum file should be written to cache"


def test_download_spectrum_data_removes_duplicates(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that duplicate rows from IAEA data are removed correctly."""
    csv_content = (
        "p_z,p_n,p_symbol,p_energy,d_z,d_n,d_symbol,bin_en,dn_de,unc_dn_de,"
        "dn_de_nu,unc_dn_de_nu,extraction_date\n"
        "38,52,Sr,0,39,51,Y,0.0,0.1,0.01,0.2,0.02,2026-04-15\n"
        "38,52,Sr,0,39,51,Y,0.0,0.1,0.01,0.2,0.02,2026-04-15\n"
        "38,52,Sr,0,39,51,Y,0.5,0.2,0.02,0.4,0.04,2026-04-15\n"
    )

    class _MockResponse:
        def read(self) -> bytes:
            return csv_content.encode("utf-8")

    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(
        "urllib.request.urlopen",
        lambda request, timeout=5: _MockResponse(),
    )

    filepath = Path(_download_spectrum_data("Sr90"))
    cached_data = np.genfromtxt(filepath, delimiter=",", skip_header=1, dtype=str)

    assert cached_data.shape[0] == 2, "Caching file should remove duplicate rows"


def test_download_spectrum_data_does_not_overwrite_existing_cache(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that an existing cache file is not overwritten by a new download."""
    cache_file = tmp_path / "90sr.csv"
    sentinel_content = "already,cached\n1,2\n"
    cache_file.write_text(sentinel_content, encoding="utf-8")

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

    filepath = Path(_download_spectrum_data("Sr90"))

    assert filepath == cache_file
    assert filepath.read_text(encoding="utf-8") == sentinel_content


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
        _download_spectrum_data("Ru106")


def test_download_spectrum_data_reports_http_errors(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test Cloudflare 403 responses are reported with a clear message."""
    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))

    def _raise_error(request: object, timeout: int = 5) -> object:
        del request, timeout
        raise urllib.error.HTTPError(
            url="https://nds.iaea.org/relnsd/v1/data",
            code=403,
            msg="Forbidden",
            hdrs=None,  # type: ignore
            fp=BytesIO(),
        )

    monkeypatch.setattr("urllib.request.urlopen", _raise_error)

    with pytest.raises(RuntimeError, match=r"HTTP error"):
        _download_spectrum_data("Ru106")


def test_download_spectrum_data_reports_cloudflare_block(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test Cloudflare 403 responses are reported with a clear message."""
    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))
    cloudflare_html = b"<html><title>Just a moment...</title>cloudflare</html>"

    def _raise_cloudflare(request: object, timeout: int = 5) -> object:
        del request, timeout
        raise urllib.error.HTTPError(
            url="https://nds.iaea.org/relnsd/v1/data",
            code=403,
            msg="Forbidden",
            hdrs=None,  # type: ignore
            fp=BytesIO(cloudflare_html),
        )

    monkeypatch.setattr("urllib.request.urlopen", _raise_cloudflare)

    with pytest.raises(RuntimeError, match=r"blocked by Cloudflare"):
        _download_spectrum_data("Ru106")


def test_load_spectrum_file_missing_cache_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that loading a missing cached spectrum raises ValueError."""
    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))

    with pytest.raises(ValueError, match=r"not found in cache"):
        _load_spectrum_file("Sr90")


def test_load_spectrum_file_contents(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test _load_spectrum_file returns the expected data."""
    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))
    cache_file = tmp_path / "90sr.csv"
    cache_file.write_text(
        "p_z,p_n,p_symbol,p_energy,d_z,d_n,d_symbol,bin_en,dn_de,unc_dn_de,"
        "dn_de_nu,unc_dn_de_nu,extraction_date\n"
        "38,52,Sr,0,39,51,Y,0.5,0.1,0.01,1.5,0.15,2026-04-15\n"
        "38,52,Sr,1,39,51,Y,1.5,0.2,0.02,2.5,0.25,2026-04-15\n"
        "38,52,Sr,0,39,51,Y,2.5,0.3,0.03,3.5,0.35,2026-04-15\n",
        encoding="utf-8",
    )

    spectrum = _load_spectrum_file("Sr90")

    assert spectrum.shape == (2, 3)
    np.testing.assert_allclose(
        spectrum,
        np.array([[0.5, 1.5, 0.15], [2.5, 3.5, 0.35]]),
    )


def test_get_antineutrino_spectrum() -> None:
    """Test loading spectrum data for a valid isotope."""
    spectrum = get_antineutrino_spectrum("Sr90")
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


def test_get_antineutrino_spectrum_invalid() -> None:
    """Test that loading invalid isotope raises ValueError."""
    with pytest.raises(ValueError, match=r"Isotope format not recognized"):
        get_antineutrino_spectrum("InvalidIsotope")


def test_get_antineutrino_spectrum_downloads_when_cache_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that get_antineutrino_spectrum triggers a download when cache is missing."""
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

    monkeypatch.setattr(
        "snf_simulations.data.iaea._download_spectrum_data", _fake_download
    )
    monkeypatch.setattr(
        "snf_simulations.data.iaea._copy_packaged_spectrum_to_cache",
        lambda _: False,
    )

    spectrum = get_antineutrino_spectrum("Sr90")

    assert called == ["90sr"]
    np.testing.assert_allclose(spectrum, np.array([[0.5, 1.5, 0.15]]))


def test_get_antineutrino_spectrum_uses_cached_file(
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
        "snf_simulations.data.iaea._download_spectrum_data",
        lambda nuclide: pytest.fail(f"Unexpected download for {nuclide}"),
    )

    spectrum = get_antineutrino_spectrum("Sr90")

    np.testing.assert_allclose(spectrum, np.array([[1.5, 2.5, 0.25]]))


def test_get_antineutrino_spectrum_copies_packaged_file(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test that packaged data is copied to cache before any download attempt."""
    monkeypatch.setenv("SNF_SIMULATIONS_CACHE_DIR", str(tmp_path))

    def _fake_copy_packaged(isotope_name: str) -> bool:
        cache_file = tmp_path / "90sr.csv"
        cache_file.write_text(
            "p_z,p_n,p_symbol,p_energy,d_z,d_n,d_symbol,bin_en,dn_de,unc_dn_de,"
            "dn_de_nu,unc_dn_de_nu,extraction_date\n"
            "38,52,Sr,0,39,51,Y,3.0,0.1,0.01,4.0,0.4,2026-04-15\n",
            encoding="utf-8",
        )
        return True

    monkeypatch.setattr(
        "snf_simulations.data.iaea._copy_packaged_spectrum_to_cache",
        _fake_copy_packaged,
    )
    monkeypatch.setattr(
        "snf_simulations.data.iaea._download_spectrum_data",
        lambda isotope_name: pytest.fail(f"Unexpected download for {isotope_name}"),
    )

    spectrum = get_antineutrino_spectrum("Sr90")

    np.testing.assert_allclose(spectrum, np.array([[3.0, 4.0, 0.4]]))
    assert (tmp_path / "90sr.csv").is_file()
