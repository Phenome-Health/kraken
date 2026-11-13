"""Moving single ids between clusters, after clustering.

A source can place an id with exactly one other id (see ``umls_minted_codes`` and ``icd10_twins``). Weighed as match
evidence, that statement changes how the ids around it cluster. Applied afterwards, it moves that id only.
"""

from __future__ import annotations

from collections.abc import Iterable


class FinalClusters:
    """Edits ``curie_to_cluster`` in place. Only the clusters of ``ids`` are looked at, so only they can be asked
    about or moved. A move may open a cluster, so ``next_cluster_id`` is the next free one."""

    def __init__(self, curie_to_cluster: dict[str, int], ids: Iterable[str], next_cluster_id: int):
        self._cluster_of = curie_to_cluster
        self.next_cluster_id = next_cluster_id
        self._members: dict[int, list[str]] = {
            curie_to_cluster[curie]: [] for curie in ids if curie in curie_to_cluster
        }
        for curie, cluster in curie_to_cluster.items():
            if cluster in self._members:
                self._members[cluster].append(curie)

    def is_alone(self, curie: str) -> bool:
        cluster = self._cluster_of.get(curie)
        return cluster is None or len(self._members[cluster]) == 1

    def members_with(self, curie: str) -> list[str]:
        """The ids of ``curie``'s cluster. An id whose evidence never reached tau has no cluster: itself alone."""
        cluster = self._cluster_of.get(curie)
        return [curie] if cluster is None else self._members[cluster]

    def join(self, curie: str, target: str) -> None:
        """Move ``curie`` out of its cluster and into ``target``'s."""
        cluster = self._cluster_of.get(target)
        if cluster is None:
            cluster, self.next_cluster_id = self.next_cluster_id, self.next_cluster_id + 1
            self._cluster_of[target] = cluster
            self._members[cluster] = [target]
        left = self._cluster_of.get(curie)
        if left is not None:
            self._members[left].remove(curie)
        self._cluster_of[curie] = cluster
        self._members[cluster].append(curie)
