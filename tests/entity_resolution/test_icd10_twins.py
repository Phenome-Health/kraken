"""A WHO ICD-10 code joins the node of the ICD-10-CM code it is after clustering, and moves nothing else (see
icd10_twins). What it does through a whole build is in tests/test_icd10_harmonizer.py."""

from kraken.entity_resolution.final_clusters import FinalClusters
from kraken.entity_resolution.guardrails import GuardrailConfig, NodeInfo
from kraken.entity_resolution.icd10_twins import attach_twins

DISEASE = NodeInfo(curie="", branches=frozenset({"disease_pheno"}))


def _attach(
    twin_of: dict[str, str],
    curie_to_cluster: dict[str, int],
    config: GuardrailConfig = GuardrailConfig(),
    facts: dict[str, NodeInfo] | None = None,
) -> set[str]:
    """The codes it moved. Every id is a disease with no taxon, but those in ``facts``."""
    moved, _next_cluster_id = attach_twins(
        twin_of,
        curie_to_cluster,
        info_of=lambda curie: (facts or {}).get(curie, DISEASE),
        guardrail_config=config,
        next_cluster_id=max(curie_to_cluster.values()) + 1,
    )
    return moved


def test_a_code_another_source_has_on_a_node_leaves_it_alone():
    """kg2 and ROBOKOP list ICD10:S12.8, a fracture of the neck, on Mondo's tracheal disorder. It goes to its
    ICD-10-CM code, and the node it leaves keeps its other ids."""
    clusters = {"MONDO:0002567": 0, "DOID:3225": 0, "ICD10:S12.8": 0, "UMLS:C0478221": 1, "ICD10CM:S12.8": 1}
    _attach({"ICD10:S12.8": "ICD10CM:S12.8"}, clusters)
    assert clusters == {"MONDO:0002567": 0, "DOID:3225": 0, "ICD10:S12.8": 1, "UMLS:C0478221": 1, "ICD10CM:S12.8": 1}


def test_the_node_an_id_left_no_longer_lists_it():
    """The next code to join that node is checked against the ids it holds now."""
    clusters = {"MONDO:0002567": 0, "ICD10:S12.8": 0, "ICD10CM:S12.8": 1}
    final = FinalClusters(clusters, clusters, next_cluster_id=2)
    final.join("ICD10:S12.8", "ICD10CM:S12.8")
    assert final.members_with("MONDO:0002567") == ["MONDO:0002567"]
    assert final.members_with("ICD10CM:S12.8") == ["ICD10CM:S12.8", "ICD10:S12.8"]


def test_a_code_joins_as_the_type_of_its_twin():
    """WHO's chapter makes a Z code a fact about care. Where its ICD-10-CM code merges through Mondo it is a disease,
    and so is the code that is that code."""
    care = NodeInfo(curie="", branches=frozenset({"clinical_obs"}))
    clusters = {"MONDO:0000001": 0, "ICD10CM:Z99": 0}
    assert _attach({"ICD10:Z99": "ICD10CM:Z99"}, clusters, facts={"ICD10:Z99": care}) == {"ICD10:Z99"}
    # and it counts as one for the next code to join that node (two ICD-10-CM codes share one only if allowed to)
    clusters = {"MONDO:0000001": 0, "ICD10CM:Z98": 0, "ICD10CM:Z99": 0}
    twins = {"ICD10:Z98": "ICD10CM:Z98", "ICD10:Z99": "ICD10CM:Z99"}
    no_one_id_rule = GuardrailConfig(enforced_prefixes=frozenset())
    assert _attach(twins, clusters, no_one_id_rule, {"ICD10:Z98": care, "ICD10:Z99": care}) == set(twins)


def test_a_code_with_a_taxon_stays_where_it_is():
    """An ICD code has no species, so one that some source gives a species is not moved, whatever its twin's node
    holds. A node's taxon is its ids' own, and this one would bring a mouse to the node."""
    mouse = NodeInfo(curie="ICD10:G61.0", branches=DISEASE.branches, taxon="NCBITaxon:10090")
    clusters = {"ICD10CM:G61.0": 0, "UMLS:C0018378": 0}
    assert not _attach({"ICD10:G61.0": "ICD10CM:G61.0"}, clusters, facts={"ICD10:G61.0": mouse})


def test_the_guardrails_decide_too():
    """Typed as an id the node holds and with no taxon, a code can only break a rule on its own prefix: with one
    WHO code to a node, a node that has ICD10:G82.4 takes no second."""
    clusters = {"ICD10CM:G80.0": 0, "ICD10:G82.4": 0}
    one_who_code = GuardrailConfig(enforced_prefixes=frozenset({"ICD10"}))
    assert not _attach({"ICD10:G80.0": "ICD10CM:G80.0"}, clusters, one_who_code)
    assert "ICD10:G80.0" not in clusters
