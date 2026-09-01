"""Shared test helpers."""

import io
import tarfile

from kraken.utils.taxonomy import TaxonNormalizer


def write_test_taxdump(directory, nodes, names) -> "object":
    """Write a minimal taxdump.tar.gz (nodes.dmp + names.dmp) into `directory` and return its path.

    Args:
        nodes: (tax_id, parent_tax_id, rank) triples
        names: tax_id -> scientific name
    """
    nodes_dmp = "".join(f"{tax_id}\t|\t{parent}\t|\t{rank}\t|\t\n" for tax_id, parent, rank in nodes)
    # Include a non-scientific name class so tests cover it being skipped
    names_dmp = "".join(
        f"{tax_id}\t|\t{name}\t|\t\t|\tscientific name\t|\n{tax_id}\t|\t{name} (common)\t|\t\t|\tcommon name\t|\n"
        for tax_id, name in names.items()
    )

    taxdump_path = directory / "taxdump.tar.gz"
    with tarfile.open(taxdump_path, "w:gz") as archive:
        for member_name, content in (("nodes.dmp", nodes_dmp), ("names.dmp", names_dmp)):
            payload = content.encode()
            info = tarfile.TarInfo(member_name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
    return taxdump_path


def build_test_taxonomy(directory, nodes, names) -> TaxonNormalizer:
    """A TaxonNormalizer over a minimal taxdump written into `directory`."""
    return TaxonNormalizer(write_test_taxdump(directory, nodes, names))
