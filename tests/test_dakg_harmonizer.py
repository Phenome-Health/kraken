"""DAKG's own records name it by an infores id the registry has deprecated."""

import yaml

from kraken.config import BUILD_CONFIG_PATH
from kraken.harmonizers.dakg import DEPRECATED_SELF_ID, DAKGHarmonizer

CURRENT_ID = "infores:drugapprovals-kp"


def _source(resource_id: str, role: str) -> dict:
    return {"resource_id": resource_id, "resource_role": f"{role}_knowledge_source"}


def _parse(*sources: dict) -> tuple[str, list[str], list[str]]:
    harmonizer = object.__new__(DAKGHarmonizer)  # __init__ builds a Biolink toolkit, which this doesn't need
    harmonizer.source_infores = CURRENT_ID
    return harmonizer._parse_trapi_sources({"sources": list(sources)})


def test_flattened_provenance_uses_the_current_id_in_either_role():
    """DAKG is the aggregator on most of its edges and the primary source on 5,183."""
    assert _parse(_source(DEPRECATED_SELF_ID, "aggregator"), _source("infores:faers", "primary"))[:2] == (
        "infores:faers",
        [CURRENT_ID],
    )
    assert _parse(_source(DEPRECATED_SELF_ID, "primary"))[0] == CURRENT_ID


def test_aggregators_copies_under_the_deprecated_id_keep_being_dropped():
    """kg2 and translator-kg-open carry 162,693 DAKG edges, every one under the old id."""
    sources = yaml.safe_load(BUILD_CONFIG_PATH.read_text())["sources"]
    assert sources["dakg"]["source_id"] == CURRENT_ID
    assert DEPRECATED_SELF_ID in sources["dakg"]["former_source_ids"]
