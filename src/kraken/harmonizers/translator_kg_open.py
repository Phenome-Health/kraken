from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import TRANSLATOR_SOURCE_ID


class TranslatorKGOpenHarmonizer(BaseHarmonizer):
    source_infores = TRANSLATOR_SOURCE_ID
    is_aggregator = True

    # Node property config
    equivalent_ids_prop = "equivalent_identifiers"
    taxon_props = {"taxon", "in_taxon"}  # both are used; unioned into a single top-level taxon list

    # NOTE: nodes also carry an `xref` field, but we deliberately do NOT fold it into equivalent_ids for now. Beyond
    # `equivalent_identifiers` its only content appears to be SMILES, case-variant InChIKeys (duplicates of
    # ones already present), and CHEMBL.TARGET curies (not 1:1 with genes?).
    # So `xref` is retained in attributes only. Revisit down the line.
