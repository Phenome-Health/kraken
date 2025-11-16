import json
import logging
import sys
from pathlib import Path
from typing import Dict, Any

import requests
from biomapper2.core.normalizer import Normalizer

from ..utils.constants import *
from ..utils.kg_io import save_to_jsonl
from ..utils.general import create_node, create_edge, create_edge_key


ACMG_API_URL = "https://actionability.clinicalgenome.org/ac/api/summ/brief"


def harmonize_acmg(input_file: Path, nodes_output: Path, edges_output: Path, biolink_version: str):
    """
    Harmonize ACMG/ClinicalGenome actionability data.

    This harmonizer processes gene-disease and variant-disease associations from the
    ACMG Clinical Genome Resource (ClinGen) Actionability Working Group.

    Args:
        input_file: Path to the downloaded JSON file from the ACMG API
        nodes_output: Path where the harmonized nodes.jsonl will be saved
        edges_output: Path where the harmonized edges.jsonl will be saved
        biolink_version: Version of Biolink Model to use for normalization
    """
    logging.info(f"Harmonizing ACMG: {input_file} -> {nodes_output}, {edges_output}")
    normalizer = Normalizer(biolink_version=biolink_version)

    nodes = dict()
    edges = dict()

    # Load the data (either from pre-downloaded file or fetch from API)
    if input_file.exists():
        logging.info(f"Loading ACMG data from {input_file}")
        with open(input_file, 'r') as f:
            data = json.load(f)
    else:
        logging.info(f"Fetching ACMG data from {ACMG_API_URL}")
        response = requests.get(ACMG_API_URL)
        response.raise_for_status()
        data = response.json()

        # Save the fetched data for future use
        input_file.parent.mkdir(parents=True, exist_ok=True)
        with open(input_file, 'w') as f:
            json.dump(data, f, indent=2)
        logging.info(f"Saved ACMG data to {input_file}")

    # Process each record
    for record in data:
        doc_id = record.get('docId')
        logging.info(f"On record {doc_id}")
        curation_type = record.get('curationType', 'Gene-Condition')
        modes_of_inheritance = record.get('modesOfInheritance', [])
        context = record['context']

        # Process both Adult and Pediatric contexts
        for context_type, context_data in context.items():
            logging.info(f"On context type {context_type}")

            # Process gene-disease associations
            genes = context_data.get('genes', [])
            for gene_info in genes:
                gene_symbol = gene_info['gene']
                gene_omim = gene_info['geneOmim']

                # Create gene node
                gene_node = _process_gene(gene_symbol, gene_omim, normalizer, record)
                nodes[gene_node[ID]] = gene_node

                # Process diseases associated with this gene
                diseases = gene_info.get('diseases', [])
                for disease_info in diseases:
                    # Create disease node from omim and preferredMondo identifiers
                    disease_node = _process_disease(disease_info, gene_omim, normalizer, record)
                    if disease_node:
                        nodes[disease_node[ID]] = disease_node

                        # Create gene-disease edge
                        edge = _create_gene_disease_edge(
                            gene_node[ID],
                            disease_node[ID],
                            doc_id,
                            context_type,
                            modes_of_inheritance,
                            record
                        )
                        if edge:
                            edge_key = create_edge_key(edge)
                            edges[edge_key] = edge

            # TODO: Add this in later. Get going with disease-gene associations first.
            # # Process variant-disease associations
            # variants = context_data.get('variants', [])
            # for variant_info in variants:
            #     variant_desc = variant_info.get('description')
            #     variant_type = variant_info.get('variantType')
            #
            #     if not variant_desc:
            #         continue
            #
            #     # Create variant node
            #     variant_node = _process_variant(variant_desc, variant_type, normalizer, record)
            #     if variant_node:
            #         nodes[variant_node[ID]] = variant_node
            #
            #         # Process diseases associated with this variant
            #         diseases = variant_info.get('diseases', [])
            #         for disease_info in diseases:
            #             # Create disease node from omim and preferredMondo identifiers
            #             disease_node = _process_disease(disease_info, normalizer, record)
            #             if disease_node:
            #                 nodes[disease_node[ID]] = disease_node
            #
            #                 # Create variant-disease edge
            #                 edge = _create_variant_disease_edge(
            #                     variant_node[ID],
            #                     disease_node[ID],
            #                     doc_id,
            #                     context_type,
            #                     modes_of_inheritance,
            #                     record
            #                 )
            #                 if edge:
            #                     edge_key = create_edge_key(edge)
            #                     edges[edge_key] = edge

    logging.info(f"Saving {len(nodes)} ACMG nodes and {len(edges)} edges")
    save_to_jsonl(nodes.values(), nodes_output, mode='w')
    save_to_jsonl(edges.values(), edges_output, mode='w')


def _process_disease(disease_info: Dict[str, Any], gene_omim: str, normalizer: Normalizer, record: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process a disease/condition into a harmonized node.
    Uses the omim and preferredMondo identifiers from the disease_info.
    """
    disease_label = disease_info.get('label', '')
    disease_omim = disease_info.get('omim')
    disease_mondo = disease_info.get('preferredMondo')

    # Build a dict of identifiers to normalize
    id_dict = {}
    if disease_mondo:
        id_dict['mondo'] = disease_mondo
    if disease_omim and disease_omim != gene_omim:  # Sometimes they incorrectly give the gene OMIM on the disease
        id_dict['omim'] = disease_omim

    if not id_dict:
        logging.error(f"No disease identifiers found for: {disease_label}")
        sys.exit(1)

    # Normalize disease identifiers to standard curies
    disease_curies_dict, _ = normalizer.get_curies(id_dict, stop_on_invalid_id=True)

    if disease_curies_dict:
        disease_curie = list(sorted(disease_curies_dict.keys(), reverse=True))[0]  # OMIM identifiers seem more accurate
        disease_iri = disease_curies_dict[disease_curie]
        equivalent_ids = list(disease_curies_dict.keys())
    else:
        # Fallback: if normalization fails, use the MONDO or OMIM directly
        logging.error(f"Could not normalize disease: {disease_label} with IDs {id_dict}. full disease item is: {disease_info}")
        sys.exit(1)

    node = create_node(
        curie=disease_curie,
        categories=['biolink:Disease'],
        equivalent_ids=equivalent_ids,
        provided_by=[CLINGEN_CURIE],
        name=disease_label,
        iri=disease_iri,
        synonyms=[disease_label] if disease_label else None
    )

    return node


def _process_gene(gene_symbol: str, gene_omim: str, normalizer: Normalizer,
                  record: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process a gene into a harmonized node.
    """
    # Build a dict of identifiers to normalize
    id_dict = {'omim': gene_omim}

    # Normalize to standard gene identifiers (HGNC, NCBIGene, etc.)
    gene_curies_dict, _ = normalizer.get_curies(id_dict, stop_on_invalid_id=True)

    if gene_curies_dict:
        gene_curie = list(gene_curies_dict.keys())[0]
        gene_iri = gene_curies_dict[gene_curie]
        equivalent_ids = list(gene_curies_dict.keys())
    else:
        # Fallback: use HGNC symbol as identifier
        logging.error(f"Could not normalize gene: {gene_symbol} - gene omim: {gene_omim}")
        sys.exit(1)

    node = create_node(
        curie=gene_curie,
        categories=['biolink:Gene'],
        equivalent_ids=equivalent_ids,
        provided_by=[CLINGEN_CURIE],
        name=gene_symbol,
        iri=gene_iri,
        synonyms=[gene_symbol]
    )

    return node


def _process_variant(variant_desc: str, variant_type: str, normalizer: Normalizer,
                     record: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process a variant into a harmonized node.
    """
    # Most variants won't have standard identifiers, so we'll use a local ID
    # Users can refine this later based on their needs
    variant_curie = f"ACMG:{record.get('docId', 'unknown')}_variant_{variant_desc.replace(' ', '_')}"

    attributes = {
        CLINGEN_CURIE: {
            'variant_description': variant_desc,
            'variant_type': variant_type,
        }
    }

    node = create_node(
        curie=variant_curie,
        categories=['biolink:SequenceVariant'],
        equivalent_ids=[variant_curie],
        provided_by=[CLINGEN_CURIE],
        name=variant_desc,
        iri=None,
        synonyms=[variant_desc],
        attributes=attributes
    )

    return node


def _create_gene_disease_edge(gene_id: str, disease_id: str, doc_id: str,
                                context_type: str, modes_of_inheritance: list,
                                record: Dict[str, Any]) -> Dict[str, Any]:
    """
    Create an edge representing a gene-disease association.
    """
    attributes = {
        CLINGEN_CURIE: {
            'doc_id': doc_id,
            'context': context_type,
            'modes_of_inheritance': modes_of_inheritance,
            'source_iri': f"https://actionability.clinicalgenome.org/ac/{context_type}/ui/stg2SummaryRpt?doc={doc_id}"
        }
    }

    if record.get('iri'):
        attributes[CLINGEN_CURIE]['source_iri_json'] = f"https://actionability.clinicalgenome.org/ac/{context_type}/api/sepio/doc/{doc_id}"

    edge = create_edge(
        subject_id=gene_id,
        object_id=disease_id,
        predicate='biolink:contributes_to',
        context_qualifier=context_type.lower(),
        primary_ks=CLINGEN_CURIE,
        knowledge_level='knowledge_assertion',
        agent_type='manual_agent',
        attributes=attributes
    )

    return edge


def _create_variant_disease_edge(variant_id: str, disease_id: str, doc_id: str,
                                   context_type: str, modes_of_inheritance: list,
                                   record: Dict[str, Any]) -> Dict[str, Any]:
    """
    Create an edge representing a variant-disease association.
    """
    attributes = {
        CLINGEN_CURIE: {
            'doc_id': doc_id,
            'context': context_type,
            'modes_of_inheritance': modes_of_inheritance,
        }
    }

    if record.get('iri'):
        attributes[CLINGEN_CURIE]['source_iri'] = record['iri']

    edge = create_edge(
        subject_id=variant_id,
        object_id=disease_id,
        predicate='biolink:related_to',  # User should refine this
        primary_ks=CLINGEN_CURIE,
        knowledge_level='knowledge_assertion',
        agent_type='manual_agent',
        attributes=attributes
    )

    return edge
