import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, Any

import requests
from biomapper2.core.normalizer import Normalizer
from bs4 import BeautifulSoup

from ..utils.constants import *
from ..utils.kg_io import save_to_jsonl
from ..utils.general import create_node, create_edge, create_edge_key


CLINGEN_ACMG_API_URL = "https://actionability.clinicalgenome.org/ac/api/summ/brief"
SCRAPED_CACHE = dict()


def harmonize_clingen(input_file: Path, nodes_output: Path, edges_output: Path, biolink_version: str):
    """
    Harmonize ClinicalGenome ACMG actionability data.

    This harmonizer processes gene-disease and variant-disease associations from the
    ACMG Clinical Genome Resource (ClinGen) Actionability Working Group.

    Args:
        input_file: Path to the downloaded JSON file from the ACMG API
        nodes_output: Path where the harmonized nodes.jsonl will be saved
        edges_output: Path where the harmonized edges.jsonl will be saved
        biolink_version: Version of Biolink Model to use for normalization
    """
    logging.info(f"Harmonizing ClinGen ACMG: {input_file} -> {nodes_output}, {edges_output}")
    normalizer = Normalizer(biolink_version=biolink_version)

    nodes = dict()
    edges = dict()

    # Load the data (either from pre-downloaded file or fetch from API)
    if input_file.exists():
        logging.info(f"Loading ClinGen ACMG data from {input_file}")
        with open(input_file, 'r') as f:
            data = json.load(f)
    else:
        logging.info(f"Fetching ClinGen ACMG data from {CLINGEN_ACMG_API_URL}")
        response = requests.get(CLINGEN_ACMG_API_URL)
        response.raise_for_status()
        data = response.json()

        # Save the fetched data for future use
        input_file.parent.mkdir(parents=True, exist_ok=True)
        with open(input_file, 'w') as f:
            json.dump(data, f, indent=2)
        logging.info(f"Saved ClinGen ACMG data to {input_file}")

    # Process each record
    for record in data:
        doc_id = record.get('docId')
        logging.info(f"On record {doc_id}")
        curation_type = record.get('curationType', 'Gene-Condition')
        modes_of_inheritance = record.get('modesOfInheritance', [])
        context = record['context']
        context_type = 'Adult'  # SKIP pediatric for now...

        if context_type in context.keys():
            context_data = context[context_type]

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

    logging.info(f"Saving {len(nodes)} ClinGen ACMG nodes and {len(edges)} edges")
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


def _create_gene_disease_edge(gene_id: str, disease_id: str, doc_id: str,
                                context_type: str, modes_of_inheritance: list,
                                record: Dict[str, Any]) -> Dict[str, Any]:
    """
    Create an edge representing a gene-disease association.
    """
    source_iri = f"https://actionability.clinicalgenome.org/ac/{context_type}/ui/stg2SummaryRpt?doc={doc_id}"
    if source_iri in SCRAPED_CACHE:
        scores = SCRAPED_CACHE[source_iri]
    else:
        scores = scrape_actionability_scores(source_iri)
        SCRAPED_CACHE[source_iri] = scores

    attributes = {
        CLINGEN_CURIE: {
            'doc_id': doc_id,
            'modes_of_inheritance': modes_of_inheritance,
            'scores': scores,
            'source_iri': source_iri,
            'source_iri_json': f"https://actionability.clinicalgenome.org/ac/{context_type}/api/sepio/doc/{doc_id}"
        }
    }

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


def scrape_actionability_scores(url, debug=False):
    """
    Scrape actionability scores from a ClinGen HTML page

    Args:
        url: Full URL to the ClinGen actionability page
             e.g., "https://actionability.clinicalgenome.org/ac/Adult/ui/stg2SummaryRpt?doc=AC102"
        debug: If True, print the HTML for inspection

    Returns:
        dict with scores, or None if scraping fails
    """
    start_time = time.time()

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        soup = BeautifulSoup(response.content, 'html.parser')

        scores = {
            'url': url,
            'outcome_intervention_pairs': []
        }

        # Find the "Final Consensus Scores" section
        # It's in a div with class "scrTable"
        scr_table = soup.find('div', class_='scrTable')

        if debug:
            print("="*80)
            print("SCORING TABLE:")
            print("="*80)
            if scr_table:
                print(scr_table.prettify())
            else:
                print("No scrTable found!")
            print("="*80)

        if not scr_table:
            elapsed_time = time.time() - start_time
            logger.warning(f"No scoring table found at {url} (took {elapsed_time:.2f}s)")
            return None

        # Find all data rows (rows with class "data")
        data_rows = scr_table.find_all('div', class_='data row')

        if debug:
            print(f"\nFound {len(data_rows)} data rows")

        for row in data_rows:
            # Extract each piece of scoring data
            # These divs have multiple classes, so we use lambda to check if class list contains our target
            oi_pair = row.find('div', class_=lambda x: x and 'oiPair' in x)
            severity = row.find('div', class_=lambda x: x and 'severity' in x and 'scrData' in x)
            likelihood = row.find('div', class_=lambda x: x and 'likelihood' in x and 'scrData' in x)
            # Note: there might not be a separate effectiveness div - check the HTML structure
            effectiveness_div = row.find_all('div', class_=lambda x: x and 'scrData' in x)
            # The third scrData div is typically effectiveness
            effectiveness = effectiveness_div[2] if len(effectiveness_div) > 2 else None
            noi = row.find('div', class_=lambda x: x and 'noi' in x and 'scrData' in x)
            total_score = row.find('div', class_=lambda x: x and 'totalScore' in x and 'scrData' in x)

            if debug:
                print(f"\nProcessing row:")
                print(f"  oi_pair: {oi_pair.get_text(strip=True) if oi_pair else 'None'}")
                print(f"  severity: {severity.get_text(strip=True) if severity else 'None'}")
                print(f"  likelihood: {likelihood.get_text(strip=True) if likelihood else 'None'}")
                print(f"  effectiveness: {effectiveness.get_text(strip=True) if effectiveness else 'None'}")
                print(f"  noi: {noi.get_text(strip=True) if noi else 'None'}")
                print(f"  total: {total_score.get_text(strip=True) if total_score else 'None'}")

            if oi_pair:
                pair_text = oi_pair.get_text(strip=True)

                # Split by " / " to separate outcome and intervention
                if ' / ' in pair_text:
                    parts = pair_text.split(' / ', 1)
                    outcome = parts[0].strip()
                    intervention = parts[1].strip()

                    pair = {
                        'outcome': outcome,
                        'intervention': intervention,
                        'scores': {
                            'severity': severity.get_text(strip=True) if severity else None,
                            'likelihood': likelihood.get_text(strip=True) if likelihood else None,
                            'effectiveness': effectiveness.get_text(strip=True) if effectiveness else None,
                            'nature_of_intervention': noi.get_text(strip=True) if noi else None,
                            'total': total_score.get_text(strip=True) if total_score else None
                        }
                    }
                    scores['outcome_intervention_pairs'].append(pair)

        elapsed_time = time.time() - start_time

        # If we didn't find any pairs, return None
        if not scores['outcome_intervention_pairs']:
            logging.warning(f"No scores found at {url} (took {elapsed_time:.2f}s)")
            return None

        logging.info(f"Successfully scraped {len(scores['outcome_intervention_pairs'])} pairs from {url} (took {elapsed_time:.2f}s)")
        logging.info(scores)
        return scores

    except Exception as e:
        elapsed_time = time.time() - start_time
        logging.error(f"Error scraping {url} after {elapsed_time:.2f}s: {e}")
        return None