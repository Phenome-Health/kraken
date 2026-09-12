# refmet.py
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any

from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import NODE_ATTRIBUTES, NODE_EQUIVALENT_IDS, NODE_NAME
from kraken.utils.kg_io import load_csv_to_dict_list, save_to_jsonl

# A cross-reference that more than this many RefMet entries share names a CLASS, not a compound. RefMet maps
# individual lipid species to the KEGG entry for their whole class -- 843 glucosylceramide species all carry
# KEGG.COMPOUND:C01190 ("Glucosylceramide"), 704 phosphatidylcholines carry C00157 -- and listing that as an
# equivalent id asserts every one of those species IS the class. In 2.1.1 that fused them: one node named
# "GlcCer 18:1;O2/13:0" ended up holding 1,262 ids. Against the 2026-08 release, "more than 5" strips exactly
# the 31 class-level KEGG ids and nothing else: 97.5% of RefMet's KEGG ids map to one entry, and no other
# prefix is ever shared by more than 5 (those small overlaps are RefMet duplicates and stereo variants, left
# for entity resolution to sort out).
MAX_ENTRIES_PER_XREF = 5
# Where a stripped class-level xref is kept instead, so the class membership isn't lost -- only the claim that
# the species is equivalent to it.
CLASS_XREFS_ATTRIBUTE = "class_level_xrefs"
MAX_CLASS_XREF_EXAMPLES = 5


class RefMetHarmonizer(BaseHarmonizer):
    """Harmonizer for RefMet CSV files"""

    attribute_props = {"super_class", "main_class", "sub_class"}
    equiv_id_props = {"pubchem_cid", "chebi_id", "hmdb_id", "lipidmaps_id", "kegg_id", "inchi_key"}

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

        for row in load_csv_to_dict_list(input_file):
            node = self._harmonize_row(row)
            if node:
                nodes[node["id"]] = node

        self._strip_class_level_xrefs(nodes)

        logging.info(f"Saving {len(nodes)} RefMet nodes")
        save_to_jsonl(nodes.values(), nodes_output, mode="w")
        save_to_jsonl([], edges_output, mode="w")  # Empty edges file

        logging.info(f"{self.source_name} harmonization complete: {len(nodes)} nodes, 0 edges")

    def _strip_class_level_xrefs(self, nodes: dict[str, dict[str, Any]]) -> None:
        """Move any xref shared by more than MAX_ENTRIES_PER_XREF entries out of equivalent_ids, in place.

        Counted over the normalized curies, after every row is harmonized, so differently spelled source
        values for the same id count together. The stripped ids are kept in an attribute (see
        CLASS_XREFS_ATTRIBUTE). A node's own RefMet id is never touched.
        """
        holders: dict[str, list[str]] = defaultdict(list)
        for node_id, node in nodes.items():
            for equiv_id in node[NODE_EQUIVALENT_IDS]:
                if equiv_id != node_id:
                    holders[equiv_id].append(node_id)
        class_xrefs = {xref: held_by for xref, held_by in holders.items() if len(held_by) > MAX_ENTRIES_PER_XREF}
        if not class_xrefs:
            return

        for xref, held_by in class_xrefs.items():
            for node_id in held_by:
                node = nodes[node_id]
                node[NODE_EQUIVALENT_IDS] = [e for e in node[NODE_EQUIVALENT_IDS] if e != xref]
                attributes = node.setdefault(NODE_ATTRIBUTES, {}).setdefault(self.source_infores, {})
                attributes.setdefault(CLASS_XREFS_ATTRIBUTE, []).append(xref)

        affected = len({node_id for held_by in class_xrefs.values() for node_id in held_by})
        ranked = sorted(class_xrefs.items(), key=lambda kv: -len(kv[1]))
        examples = [f"{xref} (x{len(held_by)}, e.g. {nodes[held_by[0]].get(NODE_NAME)!r})" for xref, held_by in ranked]
        logging.info(
            f"Stripped {len(class_xrefs)} class-level xrefs from the equivalent_ids of {affected} RefMet entries "
            f"(each was shared by more than {MAX_ENTRIES_PER_XREF} entries, so names a class, not a compound; kept "
            f"in the '{CLASS_XREFS_ATTRIBUTE}' attribute). Largest: {examples[:MAX_CLASS_XREF_EXAMPLES]}"
        )

    def _harmonize_row(self, row: dict[str, Any]) -> dict[str, Any] | None:
        # Transform the 'canonical' ID into standard curie form
        # Note: original has ' refmet_id' with leading space
        rm_curie_dict, _, _ = self.normalizer.get_curies(
            {"refmet": row[" refmet_id"]}, stop_on_invalid_id=True, fuzzy_match_vocab=False
        )
        rm_curie, rm_iri = next(iter(rm_curie_dict.items()))

        # Grab all xrefs and transform into standardized curies
        equivalent_ids = {rm_curie}
        equiv_curies_dict, _, _ = self.normalizer.get_curies(
            {prop: row[prop] for prop in self.equiv_id_props if row.get(prop)},
            stop_on_invalid_id=False,
            fuzzy_match_vocab=False,
        )
        if equiv_curies_dict:
            equivalent_ids |= set(equiv_curies_dict)

        # Put together our node
        name = row["refmet_name"]

        return self.create_node(
            curie=rm_curie,
            categories=["biolink:SmallMolecule"],
            equivalent_ids=list(equivalent_ids),
            provided_by=self.source_infores,
            name=name,
            urls=rm_iri,
            chemical_formula=row["formula"],
            exact_mass=row["exactmass"],
            attributes={k: row[k] for k in self.attribute_props},
        )
