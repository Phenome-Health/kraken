import logging
import re
import tempfile
from pathlib import Path
from typing import Any

import jsonlines

from kraken.harmonizers.base import BaseHarmonizer
from kraken.harmonizers.helpers.name_overrides import ORIGINAL_NAME_ATTRIBUTE
from kraken.utils.constants import NODE_EQUIVALENT_IDS, NODE_ID, SMILES_PREFIX
from kraken.utils.general import to_list

# ROBOKOP lists a variant's dbSNP rsid as an equivalent id of its ClinGen allele (CAID) node. They are not the
# same thing: an rsid names a POSITION, and each CAID is one ALLELE at it -- rs7944541 is carried by three
# different CAIDs. Treated as an equivalence, the rsid pulls every allele at its position into one cluster
# (92,318 rsids are shared by 2+ CAIDs, ~236k allele nodes), and even an rsid listed by a single CAID only
# claims "this position is that allele". So the pairing is emitted as an edge instead: allele member_of
# position. An HGVS expression, by contrast, spells out the change, so it IS the allele and stays an alias.
ALLELE_PREFIX = "CAID"
POSITION_PREFIX = "DBSNP"
ALLELE_TO_POSITION_PREDICATE = "biolink:member_of"  # a RefSNP is the collection of alleles at a position
POSITION_CATEGORY = "biolink:SequenceVariant"  # what Biolink uses for DBSNP ids (see its id_prefixes)

# A few thousand chemical nodes carry their structure as a bare `smiles` string (one per node). It's folded into
# equivalent_ids as a SMILES curie, which biomapper2 canonicalizes like any other -- so the same structure from
# lipidmaps or translator gets the same id.
SMILES_PROP = "smiles"

# ROBOKOP names a CAID allele after the rsid of the POSITION it sits at, so all three alleles at rs7944541 are
# called "rs7944541" -- as is the position itself. The allele's own identity is the change it makes, which the
# node carries in its HGVS expressions ("HGVS:NC_000011.8:g.30011186G>A"), so the change is appended to the
# name: "rs7944541 G>A". Append-only, like the name overrides in helpers/name_overrides.py, so the rsid is still
# a substring and text search on it still finds the allele -- but the alleles are now distinguishable from each
# other and from the position node, which keeps the bare rsid. Where the change can't be read (no HGVS, or the
# assemblies disagree on it) the allele is marked as one anyway, so it is never confusable with the position.
HGVS_PROP = "hgvs"
RSID_NAME = re.compile(r"^rs\d+$", re.IGNORECASE)
# The change at the end of an HGVS expression: a substitution, deletion, duplication, insertion or inversion,
# after the coordinate ("g.30011186G>A" -> "G>A", "g.123_124insAT" -> "insAT").
HGVS_CHANGE = re.compile(
    r"[gcnmrp]\.[\d_?+*-]+(?P<change>[ACGTUN]+>[ACGTUN]+|del[ACGTUN]*|dup[ACGTUN]*|ins[ACGTUN]+|inv[ACGTUN]*)",
    re.IGNORECASE,
)
UNKNOWN_CHANGE_SUFFIX = "allele"


class RobokopHarmonizer(BaseHarmonizer):
    is_aggregator = True

    # Node property config
    category_prop = "category"
    equivalent_ids_prop = "equivalent_identifiers"
    synonyms_props = set()
    url_prop = "url"

    # Edge property config
    publications_info_prop = "sentences"
    # Not a source we ingest, so it stays a manual exclusion. HUGE (60m edges) and we get it from
    # Translator KG anyway -- which means Translator must actually be producing edges for this to
    # be a trade rather than a loss (it silently was not, in 2.1.1).
    source_exclusions = {"infores:ubergraph"}

    def __init__(self, biolink_client, source_id: str, **kwargs):
        super().__init__(biolink_client, source_id, **kwargs)
        self.alleles_named_by_change = 0
        self.alleles_named_without_change = 0

    def harmonize(self, nodes_output: Path, edges_output: Path, **inputs: Any):
        """The default split-file harmonization, plus allele -> position edges and a node per position.

        Pairs are spooled to disk while nodes stream (ROBOKOP carries ~5M of them) and appended once both
        files are written: position nodes after the nodes, member_of edges after the edges.
        """
        with tempfile.NamedTemporaryFile("w+", suffix=".tsv", delete=True) as spool:
            self._allele_position_spool = spool
            self._allele_position_count = 0
            super().harmonize(nodes_output, edges_output, **inputs)
            self._append_positions_and_membership(nodes_output, edges_output)

    def _harmonize_node(self, node: dict[str, Any]) -> dict[str, Any]:
        node = self._name_allele_by_its_change(node)
        smiles = node.get(SMILES_PROP)
        if isinstance(smiles, str) and smiles.strip():
            equivalent_ids = to_list(node.get(self.equivalent_ids_prop))
            node = {**node, self.equivalent_ids_prop: [*equivalent_ids, f"{SMILES_PREFIX}:{smiles.strip()}"]}
        harmonized = super()._harmonize_node(node)
        allele = harmonized[NODE_ID]
        if not allele.startswith(f"{ALLELE_PREFIX}:"):
            return harmonized
        equivalent_ids = harmonized[NODE_EQUIVALENT_IDS]
        positions = [e for e in equivalent_ids if e.startswith(f"{POSITION_PREFIX}:")]
        if positions:
            harmonized[NODE_EQUIVALENT_IDS] = [e for e in equivalent_ids if not e.startswith(f"{POSITION_PREFIX}:")]
            spool = getattr(self, "_allele_position_spool", None)  # only set inside harmonize()
            for position in positions:
                if spool is not None:
                    spool.write(f"{allele}\t{position}\n")
                    self._allele_position_count += 1
        return harmonized

    def _name_allele_by_its_change(self, node: dict[str, Any]) -> dict[str, Any]:
        """Append a CAID allele's change to its rsid name (see HGVS_CHANGE). Returns the node unchanged unless it
        is an allele named after a bare rsid; the name it came with is kept in the ``original_name`` attribute."""
        curie = node.get(self.id_prop, "")
        name = node.get(self.name_prop)
        if not curie.startswith(f"{ALLELE_PREFIX}:") or not name or not RSID_NAME.match(name.strip()):
            return node
        changes = {
            match.group("change").upper()
            for expression in to_list(node.get(HGVS_PROP))
            if isinstance(expression, str) and (match := HGVS_CHANGE.search(expression))
        }
        change = changes.pop() if len(changes) == 1 else UNKNOWN_CHANGE_SUFFIX
        self.alleles_named_by_change += change != UNKNOWN_CHANGE_SUFFIX
        self.alleles_named_without_change += change == UNKNOWN_CHANGE_SUFFIX
        return {**node, self.name_prop: f"{name.strip()} {change}", ORIGINAL_NAME_ATTRIBUTE: name}

    def _append_positions_and_membership(self, nodes_output: Path, edges_output: Path) -> None:
        if not self._allele_position_count:
            return
        spool = self._allele_position_spool
        spool.flush()
        spool.seek(0)
        seen_positions: set[str] = set()
        with jsonlines.open(nodes_output, "a") as nodes, jsonlines.open(edges_output, "a") as edges:
            for line in spool:
                allele, position = line.rstrip("\n").split("\t")
                if position not in seen_positions:
                    seen_positions.add(position)
                    nodes.write(
                        self.create_node(
                            curie=position,
                            categories=[POSITION_CATEGORY],
                            provided_by=self.source_infores,
                            equivalent_ids=[position],
                            name=position.split(":", 1)[1],
                        )
                    )
                edges.write(
                    self.create_edge(
                        subject_id=allele,
                        object_id=position,
                        predicate=ALLELE_TO_POSITION_PREDICATE,
                        # ROBOKOP asserted the pairing (as a node equivalence); how it was derived isn't recorded.
                        primary_ks=self.source_infores,
                        knowledge_level="knowledge_assertion",
                        agent_type="not_provided",
                    )
                )
        logging.info(
            f"{self.source_name}: named {self.alleles_named_by_change} CAID alleles by the change they make "
            f"(e.g. 'rs7944541 G>A'), and {self.alleles_named_without_change} as '{UNKNOWN_CHANGE_SUFFIX}' where "
            f"their HGVS expressions gave no single change; the rsid alone names the POSITION they sit at"
        )
        logging.info(
            f"{self.source_name}: moved {self._allele_position_count} dbSNP rsids out of CAID allele nodes' "
            f"equivalent_ids into '{ALLELE_TO_POSITION_PREDICATE}' edges, with {len(seen_positions)} position "
            f"nodes (an rsid names a position, a CAID one allele at it, so they are not equivalent)"
        )
