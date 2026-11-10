"""The PGS scoring-file download: one network failure must not stop a build, but a lasting one must."""

import pytest
import requests

from kraken.harmonizers import pgs_catalog
from kraken.harmonizers.pgs_catalog import DOWNLOAD_ATTEMPTS, PGSCatalogHarmonizer

URL = "https://example.org/PGS1.txt.gz"


class _Response:
    def __init__(self, status_code: int = 200, chunks: tuple[bytes, ...] = (b"data",), body_error=None):
        self.status_code = status_code
        self._chunks = chunks
        self._body_error = body_error

    def __enter__(self):
        return self

    def __exit__(self, *_exc):
        return False

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code}", response=self)

    def iter_content(self, chunk_size):
        yield from self._chunks
        if self._body_error:
            raise self._body_error


@pytest.fixture
def replies(monkeypatch):
    """Scripted outcomes of requests.get, in order: a response to return, or an exception to raise."""
    outcomes: list = []

    def fake_get(_url, **_kwargs):
        outcome = outcomes.pop(0)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    monkeypatch.setattr(pgs_catalog.requests, "get", fake_get)
    monkeypatch.setattr(pgs_catalog.time, "sleep", lambda _seconds: None)
    return outcomes


def test_a_connect_timeout_is_retried(tmp_path, replies):
    """The failure that stopped the 2.4.0 build after 2 h 23 min."""
    replies += [requests.exceptions.ConnectTimeout("timed out"), _Response()]
    dest = tmp_path / "PGS1.txt.gz"
    assert PGSCatalogHarmonizer._download(URL, dest) == dest
    assert dest.read_bytes() == b"data"


def test_a_transfer_that_dies_midway_restarts_and_leaves_no_partial_file(tmp_path, replies):
    dies = _Response(chunks=(b"par",), body_error=requests.exceptions.ChunkedEncodingError("cut off"))
    replies += [dies, _Response(chunks=(b"whole",))]
    dest = tmp_path / "PGS1.txt.gz"
    PGSCatalogHarmonizer._download(URL, dest)
    assert dest.read_bytes() == b"whole"
    assert [p.name for p in tmp_path.iterdir()] == [dest.name]


def test_failure_stays_loud_when_retries_run_out(tmp_path, replies):
    replies += [requests.exceptions.ConnectTimeout("timed out")] * DOWNLOAD_ATTEMPTS
    dest = tmp_path / "PGS1.txt.gz"
    with pytest.raises(requests.exceptions.ConnectTimeout):
        PGSCatalogHarmonizer._download(URL, dest)
    assert not dest.exists()
    assert not replies  # every attempt was used


def test_a_client_error_is_not_retried(tmp_path, replies):
    replies.append(_Response(status_code=403))
    with pytest.raises(requests.HTTPError):
        PGSCatalogHarmonizer._download(URL, tmp_path / "PGS1.txt.gz")


def test_a_rate_limit_is_waited_out(tmp_path, replies):
    """429 is a client error that clears with time, so it is retried like a server error."""
    replies += [_Response(status_code=429), _Response()]
    dest = tmp_path / "PGS1.txt.gz"
    assert PGSCatalogHarmonizer._download(URL, dest) == dest


def test_404_still_means_no_scoring_file(tmp_path, replies):
    """A score without a harmonized file is normal: it just gets no gene and variant edges."""
    replies.append(_Response(status_code=404))
    assert PGSCatalogHarmonizer._download(URL, tmp_path / "PGS1.txt.gz") is None
