import re
from functools import cached_property
from typing import Any

from kraken.harmonizers.base import BaseHarmonizer
from kraken.harmonizers.helpers.source_corrections import ORIGINAL_STATEMENT_CORRECTIONS
from kraken.utils.constants import (
    EDGE_ATTRIBUTES,
    KG2_PRE_IDS_ATTR,
    KG2_STATEMENT_SEPARATOR,
    OBJ_ASPECT_QUALIFIER,
    OBJ_DIRECTION_QUALIFIER,
)

# KG2 adds each UMLS semantic type to a node's description as its own "; "-joined part, "UMLS Semantic Type: STY:T047".
# For about 2.4M records that is the whole description. They are a classification, not a description (#17).
UMLS_SEMANTIC_TYPE = re.compile(r"UMLS Semantic Type: (STY:T\d{3})")
UMLS_SEMANTIC_TYPES_ATTRIBUTE = "umls_semantic_types"
# An HTML entity's ";" is not a separator: KG2 writes "&apos;proximal tibia&apos;" (unescaped later, in create_node).
HTML_ENTITY = re.compile(r"&(#\d+|#x[0-9A-Fa-f]+|[A-Za-z][A-Za-z0-9]*);")
ENTITY_END = "\x00"


def split_umls_semantic_types(description: str) -> tuple[str, list[str]]:
    """The description without its semantic type parts, and the types. Empty parts go only where they touch a removed
    type or the ends; a run of empty parts between two pieces of text is the source's own and stays."""
    if "UMLS Semantic Type" not in description:
        return description, []
    parts = HTML_ENTITY.sub(rf"&\1{ENTITY_END}", description).split(";")
    types = [UMLS_SEMANTIC_TYPE.fullmatch(part.strip()) for part in parts]
    if not any(types):
        return description, []
    kept: list[str] = []
    run: list[str] = []  # the empty and type parts since the last piece of text
    run_has_type = False
    for part, semantic_type in zip(parts, types):
        if semantic_type or not part.strip():
            run.append(part)
            run_has_type |= bool(semantic_type)
            continue
        if kept and not run_has_type:
            kept.extend(run)
        kept.append(part)
        run, run_has_type = [], False
    text = ";".join(kept).strip().replace(ENTITY_END, ";")
    return text, sorted({semantic_type.group(1) for semantic_type in types if semantic_type})


class KG2Harmonizer(BaseHarmonizer):
    is_aggregator = True

    # Node property config
    category_prop = "all_categories"
    equivalent_ids_prop = "equivalent_curies"
    synonyms_props = {"all_names"}
    url_prop = "iri"
    rename_node_attrs = {"category": "canonical_category"}
    ignore_node_props = {"publications"}  # A lot are not very relevant, exclude

    # Edge property config
    ignore_edge_props = {"domain_range_exclusion"}
    rename_edge_attrs_or_quals = {
        "id": "kg2c_ids",
        "kg2_ids": KG2_PRE_IDS_ATTR,
        "qualified_object_direction": OBJ_DIRECTION_QUALIFIER,
        "qualified_object_aspect": OBJ_ASPECT_QUALIFIER,
    }
    # KG2 uses an invalid predicate for NCIT 'regimen_has_accepted_use_for_disease' edges - remap those
    predicate_overrides = {"biolink:drug_regulatory_status_world_wide": "biolink:treats_or_applied_or_studied_to_treat"}

    def _harmonize_node(self, node: dict[str, Any]) -> dict[str, Any]:
        # The semantic types move to an attribute. A description left empty is dropped, so the node can take a
        # lower-ranked source's real description instead.
        description = node.get(self.description_prop)
        if isinstance(description, str):  # a few are bare numbers
            text, semantic_types = split_umls_semantic_types(description)
            if semantic_types:
                node = {**node, self.description_prop: text, UMLS_SEMANTIC_TYPES_ATTRIBUTE: semantic_types}
        return super()._harmonize_node(node)

    @cached_property
    def _corrected_original_ids(self) -> frozenset[str]:
        return frozenset(curie for source, curie, _ in ORIGINAL_STATEMENT_CORRECTIONS if source == self.source_infores)

    def _harmonize_edge(self, edge: dict[str, Any]) -> dict[str, Any] | None:
        # Drop the original statements source_corrections.py names. They decide where integration puts the edge
        harmonized = super()._harmonize_edge(edge)
        attributes = harmonized.get(EDGE_ATTRIBUTES, {}).get(self.source_infores, {})
        statements = attributes.get(KG2_PRE_IDS_ATTR)
        corrected_ids = self._corrected_original_ids
        if not statements or not any(curie in statement for statement in statements for curie in corrected_ids):
            return harmonized
        kept = [statement for statement in statements if not self._is_corrected_statement(statement)]
        if not kept:
            return None
        attributes[KG2_PRE_IDS_ATTR] = kept
        return harmonized

    def _is_corrected_statement(self, statement: str) -> bool:
        parts = statement.split(KG2_STATEMENT_SEPARATOR)
        if len(parts) < 7:  # the upstream source is the 7th part
            return False
        for curie in (parts[0], parts[5]):
            key = (self.source_infores, curie, parts[-1])
            if key in ORIGINAL_STATEMENT_CORRECTIONS:
                self.applied_corrections.add(key)
                return True
        return False
