from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import MOLEPRO_INFORES, NOT_PROVIDED


class MoleProHarmonizer(BaseHarmonizer):
    source_name = "molepro"
    source_infores = MOLEPRO_INFORES
    list_delimiter = "|"
    is_aggregator = True

    # Node property config
    category_prop = "category"
    equivalent_ids_prop = "xref"
    synonyms_props = {"synonym", "trade_name", "symbol"}
    url_prop = "url"
    chemical_formula_prop = "has_chemical_formula"
    ignore_node_props = {"attributes"}
    # Fix some invalid node types that molepro uses
    category_overrides = {
        "biolink:ProteinComplex": "biolink:MacromolecularComplex",
        "biolink:MacromolecularComplexMixin": "biolink:MacromolecularComplex",
        "biolink:PairwiseGeneToGeneInteraction": "biolink:MacromolecularComplex",
        "biolink:Organism": "biolink:OrganismTaxon",
    }

    # Edge property config
    ignore_edge_props = {"attributes", "aggregator_knowledge_source"}
    predicate_overrides = {"biolink:is_active_metabolite_of": "biolink:is_metabolite_of"}
    agent_type_overrides = {"experimental_agent": {"infores:chembl": "manual_agent", "infores:bindingdb": NOT_PROVIDED}}

    # Note: Some of MolePro's names have pipes, but aren't meant to be delimited; TODO: check with them
    exclude_from_list_parsing = {"name"}
