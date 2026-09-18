"""Typed schema for build_info.json — the cross-repo KG build metadata contract.

This module is the SINGLE SOURCE OF TRUTH for the build_info.json schema. The
committed build_info.schema.json (repo root) is generated from BuildInfo by
scripts/gen_build_info_schema.py and kept in sync by tests/test_build_info.py.

Downstream consumers (Kestrel, biomapper2) keep their own model that conforms to
the same schema — see build_info.schema.json. Do not add a cross-repo import.
"""

from __future__ import annotations

from pydantic import BaseModel


class StepsRun(BaseModel):
    harmonize: bool
    integrate: bool
    postprocess: bool


class BuildInfo(BaseModel):
    # Core build provenance (always written)
    kg_version: str
    kraken_package_version: str
    biolink_version: str
    build_timestamp: str  # ISO-8601 UTC
    git_commit: str
    sources: list[str]
    steps_run: StepsRun
    build_duration_minutes: float

    # Optional analyst-facing enrichments (forward-compatible; may be absent/None)
    kg_label: str | None = None  # e.g. "kraken-full" / "kraken-lite"
    node_count: int | None = None
    edge_count: int | None = None
    # Per-source version/release ingested for this build; value is None when unknown.
    source_versions: dict[str, str | None] | None = None
