import re
from typing import Optional

from .constants import *


ILLEGAL_KEY_PATTERN = r"[^a-zA-Z0-9_\-\.:%\+\*]"


def create_edge_key(edge: dict) -> str:
    qualifiers = [
        edge.get(QUALIFIED_PREDICATE, ''),  # e.g., biolink:causes
        edge.get(QUALIFIED_DIRECTION, ''),  # e.g., increased
        edge.get(QUALIFIED_ASPECT, '')  # e.g., activity
    ]
    conglomerate_predicate = '__'.join([qualifier for qualifier in qualifiers if qualifier])
    qualifiers_str = f'{conglomerate_predicate}--' if conglomerate_predicate else ''
    subject_id = edge[SUBJECT]
    object_id = edge[OBJECT]
    predicate = edge[PREDICATE]
    primary_ks = edge[PRIMARY_KS]
    aggregator_ks = edge.get(AGGREGATOR_KS)
    aggregator_str = '--' + aggregator_ks if aggregator_ks else ''
    supporting_sources = edge.get(SUPPORTING_SOURCES)
    supporting_str = '--' + '__'.join(sorted(supporting_sources)) if supporting_sources else ''
    key_raw = f"{subject_id}--{predicate}--{qualifiers_str}{object_id}--{primary_ks}{supporting_sources}{aggregator_str}"
    return clean_key_for_arango(key_raw)


def clean_key_for_arango(key: str) -> str:
    """Remove disallowed characters to create a valid ArangoDB _key"""
    return re.sub(ILLEGAL_KEY_PATTERN, '', key)
