"""Babel reads the other sources' harmonized output, so it is harmonized last and sees every build source."""

from kraken.config import KrakenConfig
from kraken.orchestrator import harmonization_order, icd10_reads


def test_babel_is_harmonized_after_the_sources_it_reads_and_before_the_one_that_reads_it():
    assert harmonization_order({"umls", "babel", "kg2", "robokop"}) == ["kg2", "robokop", "umls", "babel"]
    assert harmonization_order({"icd10", "umls", "babel", "icd10cm"}) == ["icd10cm", "umls", "babel", "icd10"]
    assert harmonization_order({"kg2"}) == ["kg2"]


def _config(tmp_path, include=(), exclude=(), names=("babel", "kg2", "robokop", "clingen")):
    sources = {name: {"source_id": f"infores:{name}", "version": "1", "input_file": "x"} for name in names}
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


def test_icd10_reads_the_builds_other_sources_even_when_harmonized_alone(tmp_path):
    """Its harmonizer needs Babel's harmonized edges and the icd10cm source's harmonized nodes. include_sources:
    [icd10] must not hide them, or WHO's codes come out with no twin and no regard for what Babel places; an
    excluded source is not in the build."""
    names = ("babel", "icd10", "icd10cm", "kg2")
    harmonized = tmp_path / "harmonized"
    assert icd10_reads(_config(tmp_path, include=["icd10"], names=names)) == {
        "babel_edges": harmonized / "babel" / "edges.jsonl",
        "icd10cm_nodes": harmonized / "icd10cm" / "nodes.jsonl",
    }
    assert icd10_reads(_config(tmp_path, exclude=["babel"], names=names))["babel_edges"] is None
    assert icd10_reads(_config(tmp_path, exclude=["icd10cm"], names=names))["icd10cm_nodes"] is None
