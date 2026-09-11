"""Shared fixtures for entity-resolution tests."""

import pytest

from kraken.entity_resolution import build as build_mod
from kraken.entity_resolution.sri_nodenorm import NormInfo


class _OfflineNodeNorm:
    """Node normalizer stub: recognizes nothing, so the build falls back to the
    lone harmonized category (or wildcard). Keeps tests off the network. A test
    that needs specific normalizer answers monkeypatches its own stub afterward."""

    def __init__(self, *args, **kwargs):
        pass

    def resolve(self, curies, **kwargs):
        return {c: NormInfo(label=None, categories=()) for c in curies}

    def iter_cliques(self):
        return iter(())  # no normalizer cliques offline; tests needing them stub their own

    def iter_labels(self):
        return iter(())  # no per-id labels offline

    def get(self, curie):
        return None

    def close(self):
        pass


@pytest.fixture(autouse=True)
def _offline_nodenorm(monkeypatch):
    monkeypatch.setattr(build_mod, "NodeNormClient", _OfflineNodeNorm)
