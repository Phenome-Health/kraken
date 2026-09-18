"""Babel reads the other sources' harmonized output, so it is harmonized last and sees every build source."""

from kraken.config import KrakenConfig
from kraken.orchestrator import harmonization_order


def test_babel_is_harmonized_last():
    assert harmonization_order({"umls", "babel", "kg2", "robokop"}) == ["kg2", "robokop", "umls", "babel"]
    assert harmonization_order({"kg2"}) == ["kg2"]


def _config(tmp_path, include=(), exclude=()):
    sources = {
        name: {"source_id": f"infores:{name}", "version": "1", "input_file": "x"}
        for name in ("babel", "kg2", "robokop", "clingen")
    }
    return KrakenConfig(
        kraken_version="t",
        biolink_version="4.2.5",
        steps={"harmonize": True, "integrate": False, "postprocess": False},
        options={"include_sources": list(include), "exclude_sources": list(exclude)},
        base_path=str(tmp_path),
        harmonization={"output_directory": "harmonized/"},
        integration={"output_directory": "integrated/"},
        metagraph={"output_directory": "metagraphs/"},
        post_processing={"test_export": {"output_directory": "integrated/", "num_edges": 1}},
        sources=sources,
    )


def test_babel_sees_every_build_source_even_when_harmonized_alone(tmp_path):
    """A run often harmonizes one source (include_sources: [babel]) and integrates them all later, so Babel must
    look past include_sources -- but never at an excluded source, which the build won't integrate."""
    config = _config(tmp_path, include=["babel"], exclude=["clingen"])
    paths = config.harmonized_nodes_paths_of_build_sources(other_than="babel")
    assert set(paths) == {"kg2", "robokop"}
    assert all(path.name == "nodes.jsonl" and path.parent.name == name for name, path in paths.items())
