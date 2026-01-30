"""
Knowledge graph I/O utilities
"""

import csv
import gzip
import logging
import os
import shutil
import sys
import tarfile
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

import jsonlines

from kraken.utils.general import is_empty

# Increase CSV field size limit for large fields (e.g., long lists of synonyms/xrefs)
csv.field_size_limit(sys.maxsize)


def form_tarball(file_paths: list[Path], output_dir: Path, tarball_name: str | None = None):
    # Determine what to name the tarball
    if not tarball_name:
        file_name_words = [set(file_path.stem.split("_")) for file_path in file_paths]
        overlapping_words = set.intersection(*file_name_words)
        first_file_name_words = file_paths[0].stem.split("_")
        tarball_root_name = "_".join([word for word in first_file_name_words if word in overlapping_words])
        if not tarball_root_name:
            raise ValueError(
                "Could not determine name for tarball! Input files did not share any words. "
                "Please specify a 'tarball_name'."
            )
        tarball_name = f"{tarball_root_name}.tar.gz"
    else:
        if not tarball_name.endswith(".tar.gz"):
            tarball_name = f"{tarball_name}.tar.gz"

    # Actually create the tarball
    logging.info(f"Forming tarball '{tarball_name}' out of {len(file_paths)} files...")
    tarball_path = output_dir / tarball_name
    with tarfile.open(tarball_path, "w:gz") as tar:
        for file_path in file_paths:
            tar.add(file_path, arcname=file_path.name)
    logging.info(f"Done creating tarball. Saved to {tarball_path}")


def fix_repeated_prefix(type_curie: str) -> str:
    """
    Fix repeated prefixes in node/edge type CURIEs.

    Examples:
        biolink:biolink:related_to --> biolink:related_to
        biolink:biolnk:related_to --> biolink:related_to
        some_other_predicate --> some_other_predicate
    """
    parts = type_curie.split(":")
    if len(parts) <= 2:
        return type_curie
    else:
        return f"{parts[0]}:{parts[-1]}"


def _strip_key_prefixes(d: dict[str, Any]) -> dict[str, Any]:
    """Strip common prefixes like 'biolink:' from dict keys"""
    return {k.removeprefix("biolink:"): v for k, v in d.items()}


def stream_nodes_from_jsonl(nodes_file: Path) -> Iterator[dict[str, Any]]:
    """Stream nodes from JSONL file without loading into memory"""
    logging.info(f"Streaming nodes from {nodes_file}")

    with jsonlines.open(nodes_file, "r") as reader:
        for line_num, node in enumerate(reader, 1):
            if line_num % 1000000 == 0:
                logging.info(f"    at {line_num} nodes")
            yield _strip_key_prefixes(node)


def stream_edges_from_jsonl(edges_file: Path) -> Iterator[dict[str, Any]]:
    """Stream edges from JSONL file without loading into memory"""
    logging.info(f"Streaming edges from {edges_file}")

    with jsonlines.open(edges_file, "r") as reader:
        for line_num, edge in enumerate(reader, 1):
            if line_num % 5000000 == 0:
                logging.info(f"    at {line_num} edges")
            yield _strip_key_prefixes(edge)


def stream_nodes_from_tsv(
    nodes_file: Path,
    list_delimiter: str | None,
    exclude_from_list_parsing: set[str],
) -> Iterator[dict[str, Any]]:
    """Stream nodes from TSV file without loading into memory"""
    logging.info(f"Streaming nodes from {nodes_file}")

    with open(nodes_file, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for line_num, node in enumerate(reader, 1):
            if line_num % 1000000 == 0:
                logging.info(f"    at {line_num} nodes")

            # Skip blank rows
            if not any(node.values()):
                continue

            node = _strip_key_prefixes(node)

            for key, value in node.items():
                if (
                    list_delimiter
                    and not is_empty(value)
                    and key not in exclude_from_list_parsing
                    and list_delimiter in value
                ):
                    node[key] = value.split(list_delimiter)

            yield node


def stream_edges_from_tsv(
    edges_file: Path, list_delimiter: str | None, exclude_from_list_parsing: set[str]
) -> Iterator[dict[str, Any]]:
    """Stream edges from TSV file without loading into memory"""
    logging.info(f"Streaming edges from {edges_file}")

    with open(edges_file, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for line_num, edge in enumerate(reader, 1):
            if line_num % 5000000 == 0:
                logging.info(f"    at {line_num} edges")

            # Skip blank rows
            if not any(edge.values()):
                continue

            edge = _strip_key_prefixes(edge)

            for key, value in edge.items():
                if (
                    list_delimiter
                    and not is_empty(value)
                    and key not in exclude_from_list_parsing
                    and list_delimiter in value
                ):
                    edge[key] = value.split(list_delimiter)

            yield edge


def stream_mixed_jsonl(input_file: Path) -> Iterator[dict[str, Any]]:
    """Stream items from mixed JSONL file (nodes and edges together)"""
    logging.info(f"Streaming mixed JSONL from {input_file}")

    with jsonlines.open(input_file, "r") as reader:
        for line_num, item in enumerate(reader, 1):
            if line_num % 1000000 == 0:
                logging.info(f"    at {line_num} items")
            yield item


def save_to_jsonl(items: Iterable[dict], output_file_path: Path, mode: str = "w"):
    with jsonlines.open(output_file_path, mode=mode) as writer:
        writer.write_all(items)


def remove_file(file_path: Path):
    if os.path.exists(file_path):
        os.remove(file_path)


def load_equivalency_mappings(nodes_file: Path) -> dict[str, str]:
    """Load equivalency mappings for entity resolution"""
    logging.info(f"Loading equivalency mappings from {nodes_file}")

    equivalencies = {}

    for node in stream_nodes_from_jsonl(nodes_file):
        canonical_id = node["id"]
        equiv_ids = node["equivalent_ids"]

        for equiv_id in equiv_ids:
            equivalencies[equiv_id] = canonical_id

    logging.info(f"Loaded equivalencies for {len(equivalencies)} ids")
    return equivalencies


def load_csv_to_dict_list(filename: Path):
    """
    Load a CSV file into a list of dictionaries where each row becomes a dictionary
    with column headers as keys.

    Args:
        filename (str): Path to the CSV file

    Returns:
        list: List of dictionaries, one per row
    """
    records = []

    with open(filename, encoding="utf-8") as file:
        # Create a CSV reader that automatically uses the first row as headers
        csv_reader = csv.DictReader(file)

        # Convert each row to a dictionary and add to our list
        for row in csv_reader:
            records.append(dict(row))

    return records


def unzip_files(file_paths: list[str | Path | None]):
    for file_path in file_paths:
        if file_path:
            ensure_unzipped(Path(file_path))


def zip_files(file_paths: list[str | Path | None]):
    for file_path in file_paths:
        if file_path:
            ensure_gzipped(Path(file_path))


def ensure_unzipped(filepath: Path) -> Path:
    """Ensure an unzipped version of the file exists.

    Args:
        filepath: Path to the file (with or without .gz extension)

    Returns:
        Path to the unzipped file
    """
    filepath = Path(filepath)
    logging.info(f"Ensuring {filepath} is unzipped..")

    # Normalize paths for both versions
    if filepath.suffix == ".gz":
        gz_path = filepath
        unzipped_path = filepath.with_suffix("")
    else:
        unzipped_path = filepath
        gz_path = Path(str(filepath) + ".gz")

    # If unzipped exists, we're done
    if unzipped_path.exists():
        logging.info("File is already unzipped")
        return unzipped_path

    # Otherwise, need to unzip from gz
    if not gz_path.exists():
        raise FileNotFoundError(f"Neither {unzipped_path} nor {gz_path} exists")

    logging.info("File is zipped - unzipping..")
    with gzip.open(gz_path, "rb") as f_in, open(unzipped_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    return unzipped_path


def ensure_gzipped(filepath: Path, remove_original: bool = True) -> Path:
    """Ensure a gzipped version of the file exists.

    Args:
        filepath: Path to the file (with or without .gz extension)
        remove_original: If True, delete the uncompressed file after gzipping

    Returns:
        Path to the gzipped file
    """
    filepath = Path(filepath)
    logging.info(f"Ensuring {filepath} is zipped..")

    # Normalize paths for both versions
    if filepath.suffix == ".gz":
        gz_path = filepath
        unzipped_path = filepath.with_suffix("")
    else:
        unzipped_path = filepath
        gz_path = Path(str(filepath) + ".gz")

    # If gzipped exists, we're done (but maybe clean up original)
    if gz_path.exists():
        logging.info(f"File is already zipped {' - removing unzipped form' if remove_original else ''}")
        if remove_original and unzipped_path.exists():
            unzipped_path.unlink()
        return gz_path

    # Otherwise, need to zip from unzipped
    if not unzipped_path.exists():
        raise FileNotFoundError(f"Neither {unzipped_path} nor {gz_path} exists")

    logging.info("File is unzipped - zipping..")
    with open(unzipped_path, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    if remove_original:
        unzipped_path.unlink()

    return gz_path
