import html
import json
import logging
import unicodedata
from numbers import Number
from pathlib import Path
from typing import Any

import requests
import yaml

from kraken.utils.constants import (
    EDGE_OBJECT,
    EDGE_PREDICATE,
    EDGE_PRIMARY_KS,
    EDGE_QUALIFIERS,
    EDGE_SUBJECT,
    EDGE_SUPPORTING_SOURCES,
    NONE_STRINGS,
)


def create_edge_key(edge: dict) -> str:
    sep = "---"
    placeholder = "|"

    qualifiers = edge.get(EDGE_QUALIFIERS, dict())
    qualifier_strs = [f"{prop}:{qualifier}" for prop, qualifier in qualifiers.items()]
    qualifiers_str = "__".join(sorted(qualifier_strs)) if qualifier_strs else placeholder
    subject_id = edge[EDGE_SUBJECT]
    object_id = edge[EDGE_OBJECT]
    assert sep not in subject_id and sep not in object_id
    predicate = edge[EDGE_PREDICATE]
    primary_ks = edge[EDGE_PRIMARY_KS]
    supporting_sources = edge.get(EDGE_SUPPORTING_SOURCES)
    supporting_ks_str = "__".join(sorted(supporting_sources)) if supporting_sources else placeholder
    key_raw = sep.join([subject_id, predicate, object_id, qualifiers_str, primary_ks, supporting_ks_str])
    return key_raw


def clean_text(text: any) -> str:
    if not isinstance(text, str):
        # Handle weird case where some KG2 nodes have a name of True (a bool), or a description that's an int
        text = str(text)
    unescaped_text = html.unescape(text)
    cleaned_text = unicodedata.normalize("NFC", unescaped_text)
    return cleaned_text


def to_list(item: Any) -> list[Any]:
    if isinstance(item, list):
        return item
    elif isinstance(item, (set, tuple)):
        return list(item)
    elif item is None or item == "":
        return []
    else:
        return [item]


def is_empty(value: Any) -> bool:
    if isinstance(value, str) and value.lower() in NONE_STRINGS:
        return True
    elif value or isinstance(value, Number):
        return False
    else:
        return True


def load_biolink_file(url: str, biolink_version: str) -> dict:
    """Load/cache a Biolink JSON or YAML file (downloaded from a URL)"""
    project_root = Path(__file__).parents[2]
    logging.info(f"project root is: {project_root}")

    cache_dir = project_root / "cache"
    file_name = url.split("/")[-1]
    file_name_json = file_name.split(".")[0] + f"_{biolink_version}" + ".json"
    local_path = cache_dir / file_name_json
    logging.info(f"original file name is: {file_name}")
    logging.info(f"json file name is: {file_name_json}")
    logging.info(f"local path is: {local_path}")

    # Download the file if we don't already have it cached
    if not local_path.exists():
        logging.info(f"Downloading YAML file from {url}. local path is: {local_path}")
        response = requests.get(url)
        response.raise_for_status()
        if file_name.endswith(".yaml"):
            response_json = yaml.safe_load(response.text)
            print(response_json)
        else:
            response_json = response.json()

        # Cache the response
        cache_dir.mkdir(parents=True, exist_ok=True)
        with open(local_path, "w+") as cache_file:
            json.dump(response_json, cache_file, indent=2)

    # Read and return the cached JSON
    with open(local_path) as cache_file:
        contents = json.load(cache_file)
        return contents
