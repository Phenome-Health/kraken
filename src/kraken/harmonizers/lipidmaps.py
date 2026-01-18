# lipidmaps.py
import logging
from pathlib import Path
from typing import Any

from biomapper2.core.normalizer import Normalizer
from rdkit import Chem

from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.biolink_client import BiolinkClient
from kraken.utils.constants import LIPIDMAPS_ID
from kraken.utils.kg_io import save_to_jsonl


class LipidMapsHarmonizer(BaseHarmonizer):
    """Harmonizer for LIPID MAPS SDF files"""

    source_name = "lipidmaps"
    source_infores = LIPIDMAPS_ID

    attribute_props = {"CATEGORY", "MAIN_CLASS", "SUB_CLASS", "CLASS_LEVEL4", "INCHI"}
    equiv_id_props = {
        "INCHI_KEY",
        "PUBCHEM_CID",
        "CHEBI_ID",
        "KEGG_ID",
        "HMDB_ID",
        "SWISSLIPIDS_ID",
        "LIPIDBANK_ID",
        "PLANTFA_ID",
        "SMILES",
    }

    def __init__(self, biolink_client: BiolinkClient):
        super().__init__(biolink_client)
        self.normalizer = Normalizer(biolink_version=biolink_client.version)

    def harmonize(
        self,
        nodes_output: Path,
        edges_output: Path,
        *,
        input_file: Path | None = None,
        nodes_input: Path | None = None,
        edges_input: Path | None = None,
    ):
        if not input_file:
            raise ValueError(f"{self.source_name} requires input_file")

        logging.info(f"Harmonizing {self.source_name}: {input_file} -> {nodes_output}, {edges_output}")

        nodes = {}
        supplier = Chem.SDMolSupplier(str(input_file), removeHs=False)

        for molecule in supplier:
            if molecule is None:
                logging.warning("LIPID MAPS molecule failed to load")
                continue

            node = self._harmonize_molecule(molecule)
            if node:
                nodes[node["id"]] = node

        logging.info(f"Saving {len(nodes)} LIPID MAPS nodes")
        save_to_jsonl(nodes.values(), nodes_output, mode="w")
        save_to_jsonl([], edges_output, mode="w")  # Empty edges file

        logging.info(f"{self.source_name} harmonization complete: {len(nodes)} nodes, 0 edges")

    def _harmonize_molecule(self, molecule) -> dict[str, Any] | None:
        properties = molecule.GetPropsAsDict()

        # Transform the 'canonical' ID into standard curie form
        lm_curie_dict, _ = self.normalizer.get_curies({"LM_ID": properties["LM_ID"]}, stop_on_invalid_id=True)
        lm_curie, lm_iri = next(iter(lm_curie_dict.items()))

        # Grab all xrefs and transform into standardized curies
        equivalent_ids = {lm_curie}
        equiv_curies_dict, _ = self.normalizer.get_curies(
            {prop: properties[prop] for prop in self.equiv_id_props if properties.get(prop)}, stop_on_invalid_id=False
        )
        if equiv_curies_dict:
            equivalent_ids |= set(equiv_curies_dict)

        # Grab all names/synonyms
        name = properties.get("NAME")
        lm_synonyms = properties.get("SYNONYMS", "").split(";")
        other_synonyms = [name, properties.get("SYSTEMATIC_NAME"), properties.get("ABBREVIATION")]
        synonyms = {s.strip() for s in (lm_synonyms + other_synonyms) if s}

        # Collect attributes
        attributes = {k: properties[k] for k in self.attribute_props if k in properties}

        return self.create_node(
            source_infores=self.source_infores,
            curie=lm_curie,
            categories=["biolink:SmallMolecule"],
            equivalent_ids=list(equivalent_ids),
            provided_by=self.source_infores,
            name=name,
            urls=lm_iri,
            synonyms=list(synonyms),
            chemical_formula=properties.get("FORMULA"),
            exact_mass=properties.get("EXACT_MASS"),
            attributes=attributes,
        )
