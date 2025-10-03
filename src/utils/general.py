import json
import logging
import re
from pathlib import Path
from typing import Union, List, Set, Optional, Dict, Any

import requests
import yaml
from rdkit import Chem

from .constants import *


ILLEGAL_KEY_PATTERN = r"[^a-zA-Z0-9_\-\.:%\+\*]"


def create_node(curie: str,
                categories: List[str],
                equivalent_ids: List[str],
                provided_by: List[str],
                name: Optional[str] = None,
                iri: Optional[str] = None,
                synonyms: Optional[List[str]] = None,
                chemical_formula: Optional[str] = None,
                exact_mass: Optional[float] = None) -> Dict[str, Any]:
    # TODO: switch to pydantic for nodes/edges...
    assert curie and categories and equivalent_ids and provided_by

    # Assemble the node, with properties in a specific order (for convenient review)
    # TODO: add description... (check kg2 for anything else)
    node = {ID: curie}
    if name:
        node[NAME] = name
    node[CATEGORIES] = categories
    if iri:
        node[IRI] = iri
    if chemical_formula:
        node[CHEMICAL_FORMULA] = chemical_formula
    if exact_mass:
        node[EXACT_MASS] = exact_mass

    node[PROVIDED_BY] = provided_by
    node[EQUIVALENT_IDS] = equivalent_ids
    if synonyms:
        node[SYNONYMS] = synonyms

    return node


def create_edge(subject_id: str,
                object_id: str,
                predicate: str,
                primary_ks: str,
                knowledge_level: str,
                agent_type: str,
                aggregator_ks: Optional[str] = None,
                supporting_sources: Optional[List[str]] = None,
                qualified_predicate: Optional[str] = None,
                qualified_direction: Optional[str] = None,
                qualified_aspect: Optional[str] = None) -> Dict[str, Any]:
    assert subject_id and object_id and predicate and primary_ks and knowledge_level and agent_type

    # Assemble the edge, with properties in a specific order (for convenient review)
    edge = {SUBJECT: subject_id,
            OBJECT: object_id,
            PREDICATE: predicate,
            PRIMARY_KS: primary_ks,
            KNOWLEDGE_LEVEL: knowledge_level,
            AGENT_TYPE: agent_type}
    if qualified_predicate:
        edge[QUALIFIED_PREDICATE] = qualified_predicate
    if qualified_direction:
        edge[QUALIFIED_DIRECTION] = qualified_direction
    if qualified_aspect:
        edge[QUALIFIED_ASPECT] = qualified_aspect
    if supporting_sources:
        edge[SUPPORTING_SOURCES] = supporting_sources
    if aggregator_ks:
        edge[AGGREGATOR_KS] = aggregator_ks

    return edge


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
    key_raw = f"{subject_id}--{predicate}--{qualifiers_str}{object_id}--{primary_ks}{supporting_str}{aggregator_str}"
    return key_raw


def clean_key_for_arango(key: str) -> str:
    """Remove disallowed characters to create a valid ArangoDB _key"""
    return re.sub(ILLEGAL_KEY_PATTERN, '', key)


def load_biolink_file(url: str, biolink_version: str) -> dict:
    """Load/cache a Biolink JSON or YAML file (downloaded from a URL)"""
    project_root = Path(__file__).parents[2]
    logging.info(f"project root is: {project_root}")

    cache_dir = project_root / 'cache'
    file_name = url.split('/')[-1]
    file_name_json = file_name.split('.')[0] + f"_{biolink_version}" + '.json'
    local_path = cache_dir / file_name_json
    logging.info(f"original file name is: {file_name}")
    logging.info(f"json file name is: {file_name_json}")
    logging.info(f"local path is: {local_path}")

    # Download the file if we don't already have it cached
    if not local_path.exists():
        logging.info(f"Downloading YAML file from {url}. local path is: {local_path}")
        response = requests.get(url)
        response.raise_for_status()
        if file_name.endswith('.yaml'):
            response_json = yaml.safe_load(response.text)
            print(response_json)
        else:
            response_json = response.json()

        # Cache the response
        cache_dir.mkdir(parents=True, exist_ok=True)
        with open(local_path, 'w+') as cache_file:
            json.dump(response_json, cache_file, indent=2)

    # Read and return the cached JSON
    with open(local_path, 'r') as cache_file:
        contents = json.load(cache_file)
        return contents


def get_canonical_smiles(smiles_string: str) -> Optional[str]:
    """
    Convert a SMILES string to its canonical form using RDKit

    Args:
        smiles_string (str): Input SMILES string

    Returns:
        str: Canonical SMILES string, or None if invalid
    """
    try:
        # Parse the SMILES string into a molecule object
        mol = Chem.MolFromSmiles(smiles_string)

        # Check if parsing was successful
        if mol is None:
            logging.warning(f"Invalid SMILES - could not get canonical version: {smiles_string}")
            return None

        # Generate canonical SMILES
        canonical_smiles = Chem.MolToSmiles(mol)
        return canonical_smiles

    except Exception as e:
        logging.error(f"Error processing SMILES: {e}")
        return None
