from src.harmonizers.base import BaseHarmonizer
from ..utils.constants import MOLEPRO_INFORES, AGGREGATOR_KS


class MoleProHarmonizer(BaseHarmonizer):
    source_name = "molepro"
    source_infores = MOLEPRO_INFORES
    list_delimiter = "|"

    # Node property config
    category_prop = "category"
    equivalent_ids_prop = "xref"
    synonyms_props = {"synonym", "trade_name", "symbol"}
    url_prop = "url"
    chemical_formula_prop = "has_chemical_formula"
    ignore_node_props = {"attributes"}

    # Edge property config
    ignore_edge_props = {"attributes", AGGREGATOR_KS}

    # Note: Some of MolePro's names have pipes, but aren't meant to be delimited; TODO: check with them
    exclude_from_list_parsing = {"name"}
