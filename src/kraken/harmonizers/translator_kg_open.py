import json
import shutil
from pathlib import Path
from typing import Any

from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import SMILES_PREFIX
from kraken.utils.general import to_list

RELEASE_METADATA_FILENAME = "release-metadata.json"
BABEL_VERSION_FIELD = "babel_version"


def release_babel_version(nodes_path: Path) -> str:
    """The Babel release a translator-kg-open release was normalized with, from the metadata file beside its nodes:
    the release's own, or the copy harmonizing puts beside the harmonized nodes, which entity resolution reads.

    Entity resolution decides from this whether the source's lists only copy our Babel's cliques. A wrong guess either
    way moves ids (see ``babel_outliers.drop_corroborated``), so a missing file, or a field that isn't a plain
    non-blank string, is an error: compared as it stands, a null, a number or "2026jul22 " would just read as
    "another release". Any other string is taken as it is, and entity resolution logs a release that isn't ours.
    """
    path = Path(nodes_path).with_name(RELEASE_METADATA_FILENAME)
    version = json.loads(path.read_text()).get(BABEL_VERSION_FIELD)
    if not isinstance(version, str) or not version or version != version.strip():
        raise ValueError(f"{path} does not say which Babel release it used ({BABEL_VERSION_FIELD}: {version!r})")
    return version


class TranslatorKGOpenHarmonizer(BaseHarmonizer):
    is_aggregator = True

    # Node property config. translator's equivalent_identifiers are Node Normalizer cliques, some broad enough to
    # bridge distinct entities (e.g. a glycolipid class fused with its sphingomyelin member). They used to be
    # ignored for that reason; entity resolution now guards against it instead -- ids a list holds in bulk get
    # only a weak link (ERWeights.max_ids_per_prefix), all Babel-derived lists count once between them, and the
    # guardrails split what shouldn't merge -- so the mappings Babel lacks are worth having.
    equivalent_ids_prop = "equivalent_identifiers"
    taxon_props = {"taxon", "in_taxon"}  # both are used; unioned into a single top-level taxon list
    synonyms_props = {"synonym", "full_name", "symbol"}

    # Nodes also carry an `xref` field. Beyond `equivalent_identifiers` its content is SMILES, case-variant
    # InChIKeys (duplicates of ones already present) and CHEMBL.TARGET curies (not 1:1 with genes). Only the
    # SMILES are folded into equivalent_ids -- one per node, which biomapper2 canonicalizes like any other curie
    # (so the same structure from lipidmaps gets the same id). The whole `xref` is still kept in attributes.
    xref_prop = "xref"

    # involved_in doesn't exist in biolink v4.2.5 (or the more recent v4.4.4)
    predicate_overrides = {"biolink:involved_in": "biolink:actively_involved_in"}

    def harmonize(self, nodes_output: Path, edges_output: Path, **inputs: Any):
        """The default split-file harmonization, plus the release's metadata copied beside the output: entity
        resolution reads only harmonized files, and needs the Babel release from it (see ``release_babel_version``).
        Checked first, so a bad file fails before the harmonize, not after it."""
        nodes_input = Path(inputs["nodes_input"])
        release_babel_version(nodes_input)
        super().harmonize(nodes_output, edges_output, **inputs)
        shutil.copyfile(
            nodes_input.with_name(RELEASE_METADATA_FILENAME), Path(nodes_output).with_name(RELEASE_METADATA_FILENAME)
        )

    def _harmonize_node(self, node: dict[str, Any]) -> dict[str, Any]:
        smiles = [
            xref
            for xref in to_list(node.get(self.xref_prop))
            if isinstance(xref, str) and xref.upper().startswith(f"{SMILES_PREFIX}:")
        ]
        if smiles:
            node = {**node, self.equivalent_ids_prop: to_list(node.get(self.equivalent_ids_prop)) + smiles}
        return super()._harmonize_node(node)
