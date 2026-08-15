from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import TRANSLATOR_SOURCE_ID


class TranslatorKGOpenHarmonizer(BaseHarmonizer):
    source_infores = TRANSLATOR_SOURCE_ID
    is_aggregator = True

    # Node property config. We deliberately do NOT use translator's equivalent_identifiers for node merging: as
    # an SRI Node Normalizer-based aggregator, its nodes carry broad equivalence sets that bridge distinct
    # entities (e.g. the glycolipid class fused with its sphingomyelin member), over-merging the graph. Setting
    # this to "" makes each node keep only its own primary id; the raw equivalent_identifiers still land in
    # attributes for debugging. (This replaces the can_merge_existing_nodes=False workaround, which could evict
    # bridging ids from the node set and orphan edges -- to be handled properly in the ER overhaul.)
    equivalent_ids_prop = ""
    taxon_props = {"taxon", "in_taxon"}  # both are used; unioned into a single top-level taxon list

    # NOTE: nodes also carry an `xref` field, but we deliberately do NOT fold it into equivalent_ids for now. Beyond
    # `equivalent_identifiers` its only content appears to be SMILES, case-variant InChIKeys (duplicates of
    # ones already present), and CHEMBL.TARGET curies (not 1:1 with genes?).
    # So `xref` is retained in attributes only. Revisit down the line.

    # involved_in doesn't exist in biolink v4.2.5 (or the more recent v4.4.4)
    predicate_overrides = {"biolink:involved_in": "biolink:actively_involved_in"}
