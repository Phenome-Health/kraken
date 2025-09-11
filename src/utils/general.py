import json
import logging
import re
from pathlib import Path

import requests
import yaml

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
