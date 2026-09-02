import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from kraken.build_info import BuildInfo, StepsRun

SCHEMA_PATH = Path(__file__).resolve().parents[1] / "build_info.schema.json"


def test_build_info_roundtrips_required_fields():
    info = BuildInfo(
        kg_version="2026.06.0",
        kraken_package_version="0.1.0",
        biolink_version="4.2.5",
        build_timestamp="2026-06-15T00:00:00+00:00",
        git_commit="abc123",
        sources=["kg2", "umls"],
        steps_run=StepsRun(harmonize=True, integrate=True, postprocess=True),
        build_duration_minutes=42.0,
    )
    dumped = info.model_dump()
    assert dumped["kg_version"] == "2026.06.0"
    assert dumped["sources"] == ["kg2", "umls"]
    assert dumped["steps_run"]["harmonize"] is True
    # Optional analyst enrichments default to None and are present in the contract
    assert dumped["kg_label"] is None
    assert dumped["node_count"] is None


def test_build_info_rejects_missing_required_field():
    with pytest.raises(ValidationError):
        BuildInfo(kg_version="x")  # type: ignore[call-arg]


def test_committed_schema_matches_model():
    """build_info.schema.json must be regenerated when BuildInfo changes."""
    committed = json.loads(SCHEMA_PATH.read_text())
    assert (
        committed == BuildInfo.model_json_schema()
    ), "build_info.schema.json is stale — run scripts/gen_build_info_schema.py"


def test_write_build_info_produces_valid_file(tmp_path, monkeypatch):
    """_write_build_info writes a BuildInfo-valid JSON with sorted sources (no set crash)."""
    from kraken.orchestrator import KrakenBuildOrchestrator

    class _Cfg:
        kraken_version = "2026.06.0"
        biolink_version = "4.2.5"
        sources_to_use = {"umls", "kg2"}  # a SET — must be serialized as sorted list
        kg_label = "kraken-lite"
        steps = type("S", (), {"harmonize": True, "integrate": True, "postprocess": False})()
        integrated_dir = tmp_path

    orch = KrakenBuildOrchestrator.__new__(KrakenBuildOrchestrator)
    orch.config = _Cfg()  # type: ignore[attr-defined]
    orch._write_build_info(120.0)

    written = json.loads((tmp_path / "build_info.json").read_text())
    BuildInfo.model_validate(written)  # raises if invalid
    assert written["sources"] == ["kg2", "umls"]  # sorted
    assert written["kg_label"] == "kraken-lite"
    assert written["steps_run"]["postprocess"] is False
