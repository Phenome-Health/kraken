import html
import json
import logging
import unicodedata
from pathlib import Path
from typing import Union, List, Set, Optional, Dict, Any, Collection

import requests
import yaml

from .constants import *


def create_edge_key(edge: dict) -> str:
    sep = '---'
    placeholder = '|'
    qualifiers = [
        edge.get(QUALIFIED_PREDICATE, ''),  # e.g., biolink:causes
        edge.get(OBJ_DIRECTION_QUALIFIER, ''),  # e.g., increased
        edge.get(OBJ_ASPECT_QUALIFIER, ''),  # e.g., activity
        edge.get(CONTEXT_QUALIFIER, '')  # e.g., pediatric
    ]
    conglomerate_qualifiers = '__'.join([qualifier for qualifier in qualifiers if qualifier])
    qualifiers_str = conglomerate_qualifiers if conglomerate_qualifiers else placeholder
    subject_id = edge[SUBJECT]
    object_id = edge[OBJECT]
    assert sep not in subject_id and sep not in object_id
    predicate = edge[PREDICATE]
    primary_ks = edge[PRIMARY_KS]
    supporting_sources = edge.get(SUPPORTING_SOURCES)
    supporting_str = '__'.join(sorted(supporting_sources)) if supporting_sources else placeholder
    key_raw = sep.join([subject_id, predicate, object_id, qualifiers_str, primary_ks, supporting_str])
    return key_raw


def clean_text(text: any) -> str:
    if not isinstance(text, str):
        # Handle weird case where some KG2 nodes have a name of True (a bool), or a description that's an int
        text = str(text)
    unescaped_text = html.unescape(text)
    cleaned_text = unicodedata.normalize('NFC', unescaped_text)
    return cleaned_text


def to_list(item: Any) -> List[Any]:
    if isinstance(item, list):
        return item
    elif isinstance(item, (set, tuple)):
        return list(item)
    elif item is None:
        return []
    else:
        return [item]


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
