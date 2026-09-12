"""Dropping edges that belong to a source we ingest directly.

Two halves: build_config's ``drop_from_other_sources`` decides WHICH sources each harmonizer must
drop, and ``_edge_knowledge_sources`` decides whether a given edge belongs to one of them. Both
primary and aggregator provenance count -- an edge that merely passed THROUGH an excluded source is
just as second-hand as one it asserted.
"""

import pytest

from kraken.config import KrakenConfig
from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import KRAKEN_SOURCE_ID
from tests.helpers import stub_normalization


class _Harmonizer(BaseHarmonizer):
    pass


def _harmonizer(exclusions: set[str]) -> BaseHarmonizer:
    """A bare harmonizer (BaseHarmonizer.__init__ builds a Biolink toolkit, which tests don't need)."""
    instance = object.__new__(_Harmonizer)
    instance.source_exclusions = exclusions
    return instance


def _trapi(**roles: str) -> dict:
    """An edge with TRAPI-style provenance: keyword is the resource id, value is its role."""
    return {"sources": [{"resource_id": rid, "resource_role": role} for rid, role in roles.items()]}


# ------------------------------------------------------------------ reading an edge's provenance


def test_trapi_primary_and_aggregator_both_count():
    h = _harmonizer(set())
    edge = _trapi(
        **{
            "infores:multiomics-clinicaltrials": "primary_knowledge_source",
            "infores:aact": "aggregator_knowledge_source",
        }
    )
    assert h._edge_knowledge_sources(edge) == {"infores:multiomics-clinicaltrials", "infores:aact"}


def test_trapi_supporting_data_source_does_not_count():
    """A supporting_data_source contributed underlying data, not the assertion. Excluding on it
    would drop edges that source never published."""
    h = _harmonizer(set())
    edge = _trapi(**{"infores:gtex": "primary_knowledge_source", "infores:chembl": "supporting_data_source"})
    assert h._edge_knowledge_sources(edge) == {"infores:gtex"}


def test_flat_props_primary_and_aggregator_both_count():
    h = _harmonizer(set())
    edge = {
        "primary_knowledge_source": "infores:gtex",
        "aggregator_knowledge_source": ["infores:robokop-kg", "infores:other"],
    }
    assert h._edge_knowledge_sources(edge) == {"infores:gtex", "infores:robokop-kg", "infores:other"}


def test_edge_with_unreadable_provenance_yields_nothing():
    """Empty means 'unknown', and an unknown edge must be KEPT.

    The regression this guards: the old filter asked ``primary_kses.issubset(exclusions)``, and
    ``set().issubset(anything)`` is True -- so a source whose provenance it could not read had
    EVERY edge excluded. That emptied translator-kg-open (TRAPI sources, no flat primary_ks prop)
    of all 19.8M of its edges.
    """
    h = _harmonizer({"infores:multiomics-clinicaltrials"})
    for edge in ({}, {"subject": "A:1", "object": "B:2"}, {"sources": []}):
        sources = h._edge_knowledge_sources(edge)
        assert sources == set()
        assert not (h.source_exclusions & sources)  # ...so the filter keeps it


# ------------------------------------------------------------------ the filter decision


@pytest.mark.parametrize(
    "edge, excluded",
    [
        (_trapi(**{"infores:multiomics-clinicaltrials": "primary_knowledge_source"}), True),
        (_trapi(**{"infores:multiomics-clinicaltrials": "aggregator_knowledge_source"}), True),
        (_trapi(**{"infores:multiomics-clinicaltrials": "supporting_data_source"}), False),
        (_trapi(**{"infores:aact": "primary_knowledge_source"}), False),
        ({"primary_knowledge_source": "infores:multiomics-clinicaltrials"}, True),
        ({"aggregator_knowledge_source": ["infores:multiomics-clinicaltrials"]}, True),
        ({"primary_knowledge_source": "infores:gtex"}, False),
    ],
)
def test_exclusion_decision(edge, excluded):
    h = _harmonizer({"infores:multiomics-clinicaltrials"})
    assert bool(h.source_exclusions & h._edge_knowledge_sources(edge)) is excluded


# ------------------------------------------------------------------ build_config -> exclusions


def _config(**source_flags: bool) -> KrakenConfig:
    return KrakenConfig(
        biolink_version="4.2.5",
        kraken_version="0.0.0",
        harmonization={"output_directory": "h"},
        integration={"output_directory": "i"},
        metagraph={"output_directory": "m"},
        steps={"harmonize": True, "integrate": True, "postprocess": False},
        options={},
        sources={
            name: {"source_id": f"infores:{name}", "version": "1", "drop_from_other_sources": flag}
            for name, flag in source_flags.items()
        },
    )


def test_flagged_source_is_excluded_from_every_other_source():
    config = _config(ctkg=True, dakg=True, kg2=False, robokop=False)
    assert config.auto_source_exclusions("kg2") == {"infores:ctkg", "infores:dakg"}
    assert config.auto_source_exclusions("robokop") == {"infores:ctkg", "infores:dakg"}


def test_a_flagged_source_never_excludes_itself():
    """Otherwise the source we ingest directly would drop its own edges."""
    config = _config(ctkg=True, dakg=True)
    assert config.auto_source_exclusions("ctkg") == {"infores:dakg"}
    assert config.auto_source_exclusions("dakg") == {"infores:ctkg"}


def test_unflagged_sources_contribute_nothing():
    config = _config(kg2=False, robokop=False)
    assert config.auto_source_exclusions("kg2") == set()


def test_only_sources_in_this_build_count():
    """Dropping a source's second-hand copies while its own ingest is switched off would delete
    those edges from the graph entirely, so an excluded/unselected source contributes nothing."""
    config = KrakenConfig(
        biolink_version="4.2.5",
        kraken_version="0.0.0",
        harmonization={"output_directory": "h"},
        integration={"output_directory": "i"},
        metagraph={"output_directory": "m"},
        steps={"harmonize": True, "integrate": True, "postprocess": False},
        options={"exclude_sources": ["ctkg"]},
        sources={
            "ctkg": {"source_id": "infores:ctkg", "version": "1", "drop_from_other_sources": True},
            "dakg": {"source_id": "infores:dakg", "version": "1", "drop_from_other_sources": True},
            "kg2": {"source_id": "infores:kg2", "version": "1"},
        },
    )
    assert config.auto_source_exclusions("kg2") == {"infores:dakg"}  # ctkg is not in this build


def test_flag_defaults_off():
    config = _config(kg2=False)
    assert config.sources["kg2"].drop_from_other_sources is False


# ------------------------------------------------------------------ KRAKEN as aggregator of direct ingests


class _DirectSource(BaseHarmonizer):
    is_aggregator = False


class _AggregatorSource(BaseHarmonizer):
    is_aggregator = True


def _edge_from(cls, aggregator_ks=None):
    h = stub_normalization(object.__new__(cls))
    h.source_infores = "infores:test"
    h.predicate_overrides = {}
    return h.create_edge(
        subject_id="A:1",
        object_id="B:2",
        predicate="biolink:related_to",
        primary_ks="infores:test",
        knowledge_level="knowledge_assertion",
        agent_type="manual_agent",
        aggregator_ks=aggregator_ks,
    )


def test_direct_source_edges_name_kraken_as_aggregator():
    """Otherwise an edge we ingested straight from a source can't be told apart from kg2's copy of it."""
    assert _edge_from(_DirectSource)["aggregator_knowledge_source"] == [KRAKEN_SOURCE_ID]


def test_kraken_is_appended_last_to_an_existing_chain():
    """A direct source whose own records name an upstream aggregator (dakg: FAERS -> drugapprovals)
    keeps that chain, with KRAKEN as the final hop."""
    edge = _edge_from(_DirectSource, aggregator_ks=["infores:multiomics-drugapprovals"])
    assert edge["aggregator_knowledge_source"] == ["infores:multiomics-drugapprovals", KRAKEN_SOURCE_ID]


def test_kraken_is_not_duplicated():
    edge = _edge_from(_DirectSource, aggregator_ks=[KRAKEN_SOURCE_ID])
    assert edge["aggregator_knowledge_source"] == [KRAKEN_SOURCE_ID]


def test_aggregator_source_edges_do_not_get_kraken():
    """kg2/ROBOKOP/Translator record themselves; leaving KRAKEN off is what distinguishes their copies."""
    assert _edge_from(_AggregatorSource, aggregator_ks=["infores:rtx-kg2"])["aggregator_knowledge_source"] == [
        "infores:rtx-kg2"
    ]
    assert "aggregator_knowledge_source" not in _edge_from(_AggregatorSource)
