from typing import Any

from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import SMILES_PREFIX
from kraken.utils.general import to_list


class TranslatorKGOpenHarmonizer(BaseHarmonizer):
    is_aggregator = True

    # Node property config. translator's equivalent_identifiers are Node Normalizer cliques, some broad enough to
    # bridge distinct entities (e.g. a glycolipid class fused with its sphingomyelin member). They used to be
    # ignored for that reason; entity resolution now guards against it instead -- ids a list holds in bulk get
    # only a weak link (ERWeights.max_ids_per_prefix), all Babel-derived lists count once between them, and the
    # guardrails split what shouldn't merge -- so the mappings Babel lacks are worth having.
    equivalent_ids_prop = "equivalent_identifiers"
    taxon_props = {"taxon", "in_taxon"}  # both are used; unioned into a single top-level taxon list
    synonyms_props = {"synonym", "full_name", "symbol"}

    # Nodes also carry an `xref` field. Beyond `equivalent_identifiers` its content is SMILES, case-variant
    # InChIKeys (duplicates of ones already present) and CHEMBL.TARGET curies (not 1:1 with genes). Only the
    # SMILES are folded into equivalent_ids -- one per node, which biomapper2 canonicalizes like any other curie
    # (so the same structure from lipidmaps gets the same id). The whole `xref` is still kept in attributes.
    xref_prop = "xref"

    # involved_in doesn't exist in biolink v4.2.5 (or the more recent v4.4.4)
    predicate_overrides = {"biolink:involved_in": "biolink:actively_involved_in"}

    def _harmonize_node(self, node: dict[str, Any]) -> dict[str, Any]:
        smiles = [
            xref
            for xref in to_list(node.get(self.xref_prop))
            if isinstance(xref, str) and xref.upper().startswith(f"{SMILES_PREFIX}:")
        ]
        if smiles:
            node = {**node, self.equivalent_ids_prop: to_list(node.get(self.equivalent_ids_prop)) + smiles}
        return super()._harmonize_node(node)
