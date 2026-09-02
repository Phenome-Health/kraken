"""Branch-family map for the "one branch" guardrail (plan §3).

``branches(node)`` = union over the node's categories of each category's family
set (permissive within a node). A cluster is valid iff the intersection of
``branches`` across its nodes is non-empty (strict across nodes). A node typed
only ``NamedThing``/``Entity`` — or untyped — is a **wildcard**: its branch set
is "ALL", so it never causes a violation and never holds two branches together.

The map is set-valued because Biolink has multiple inheritance and mixins (e.g.
``Cell`` is both anatomy and organism; ``NucleicAcidEntity`` is both chemical and
gene_protein). It is a **curated artifact**, not derivable from the hierarchy —
see ``config/entity_resolution/branch_families.yaml``.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable
from pathlib import Path

import yaml

from kraken.utils.constants import PROJECT_ROOT

# Categories that carry no branch information: a node typed only with these (or
# untyped) is a wildcard. Category values are always biolink-prefixed after
# harmonization, so only the prefixed forms are needed here.
WILDCARD_CATEGORIES: frozenset[str] = frozenset({"biolink:NamedThing", "biolink:Entity"})

DEFAULT_FAMILIES_PATH = PROJECT_ROOT / "config" / "entity_resolution" / "branch_families.yaml"

# Sentinel meaning "compatible with every family" (wildcard nodes).
ALL_FAMILIES = frozenset({"__ALL__"})


class BranchFamilies:
    """Curated category -> family-set map, with node/cluster branch logic."""

    def __init__(self, category_families: dict[str, frozenset[str]]):
        self.category_families = category_families
        self.families: set[str] = set()
        for fams in category_families.values():
            self.families |= fams

    @classmethod
    def load(cls, path: str | Path | None = None) -> BranchFamilies:
        """Load the ``category -> [families]`` map. A category mapped to an empty
        list is a wildcard (contributes no constraint). Set-valued (a category may
        list several families) to express bridge categories."""
        path = Path(path) if path is not None else DEFAULT_FAMILIES_PATH
        data = yaml.safe_load(Path(path).read_text()) or {}
        categories = data.get("categories", {})
        category_families = {category: frozenset(fams or []) for category, fams in categories.items()}
        return cls(category_families)

    def families_for_category(self, category: str) -> frozenset[str]:
        """Family set for one category (biolink-prefixed, as harmonization emits).
        Wildcard categories -> ALL. Unknown concrete categories -> ALL (treated as
        wildcard) with a warning, so a curation gap never manufactures a spurious
        guardrail violation."""
        if category in WILDCARD_CATEGORIES:
            return ALL_FAMILIES
        fams = self.category_families.get(category)
        if fams is None:
            logging.warning("Category %r not in branch-family map; treating as wildcard.", category)
            return ALL_FAMILIES
        return fams

    def branches(self, categories: Iterable[str] | None) -> frozenset[str]:
        """branches(node) = union over its categories (permissive within a node).

        Untyped or wildcard-only nodes return ALL.
        """
        categories = [c for c in (categories or []) if c]
        if not categories:
            return ALL_FAMILIES
        result: set[str] = set()
        for category in categories:
            fams = self.families_for_category(category)
            if fams is ALL_FAMILIES or fams == ALL_FAMILIES:
                # a wildcard category doesn't narrow the union, but if EVERY
                # category is wildcard the node is a wildcard (handled below)
                continue
            result |= fams
        if not result:
            return ALL_FAMILIES
        return frozenset(result)

    @staticmethod
    def cluster_branches(node_branches: Iterable[frozenset[str]]) -> frozenset[str] | None:
        """Intersection of branches across nodes (strict across nodes).

        Wildcard nodes (ALL) don't constrain the intersection. Returns the
        surviving branch set, or ``None`` if empty (a violation). An all-wildcard
        cluster returns ALL.
        """
        result: set[str] | None = None
        for nb in node_branches:
            if nb is ALL_FAMILIES or nb == ALL_FAMILIES:
                continue
            if result is None:
                result = set(nb)
            else:
                result &= nb
            if not result:
                return None
        if result is None:
            return ALL_FAMILIES
        return frozenset(result)

    def is_valid_cluster(self, node_category_lists: Iterable[Iterable[str] | None]) -> bool:
        """A cluster satisfies the one-branch guardrail iff the intersection of
        its nodes' branch sets is non-empty."""
        branches = [self.branches(cats) for cats in node_category_lists]
        return self.cluster_branches(branches) is not None
