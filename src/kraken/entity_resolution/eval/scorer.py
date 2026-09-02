"""Pairwise precision/recall scoring of a clustering against ground truth.

Ground truth is expressed as *partitions* (see
``config/entity_resolution/ground_truth_seed.jsonl``): each case names a set of
CURIEs grouped into the true distinct entities. From those partitions we derive
constrained pairs:

* **must_link**  — two CURIEs in the *same* group of a case must land in the
  same predicted cluster.
* **cannot_link** — two CURIEs in *different* groups of the *same* case must land
  in *different* predicted clusters.

Pairs whose two CURIEs come from different cases carry no constraint (the ground
truth is local to a case). CURIEs missing from the predicted membership map are
reported as uncovered and excluded from the metrics rather than counted as
errors.

Precision/recall are computed over must_link pairs, treating "predicted same
cluster" as the positive call:

    TP = must_link pairs predicted same         FN = must_link pairs predicted apart
    FP = cannot_link pairs predicted same       TN = cannot_link pairs predicted apart
    precision = TP / (TP + FP)   recall = TP / (TP + FN)

This matches the plan's "report pairwise precision/recall as a build statistic".
"""

from __future__ import annotations

import argparse
import json
import logging
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from itertools import combinations, product
from pathlib import Path

from kraken.utils.constants import PROJECT_ROOT

Membership = Mapping[str, object]  # curie -> cluster id (any hashable)

DEFAULT_GROUND_TRUTH_PATH = PROJECT_ROOT / "config" / "entity_resolution" / "ground_truth_seed.jsonl"


@dataclass(frozen=True)
class GoldCase:
    """One ground-truth case: an id-partition into true distinct entities.

    ``groups`` is a list of frozensets; ids within a group must_link, ids across
    groups (within this case) cannot_link. Singleton groups (e.g. each member of
    a numbered name family) contribute only cannot_link constraints.
    """

    case: str
    kind: str
    groups: tuple[frozenset[str], ...]
    needs_review: frozenset[str] = frozenset()

    def must_link_pairs(self) -> Iterable[tuple[str, str]]:
        for group in self.groups:
            yield from combinations(sorted(group), 2)

    def cannot_link_pairs(self) -> Iterable[tuple[str, str]]:
        for g1, g2 in combinations(self.groups, 2):
            for a, b in product(sorted(g1), sorted(g2)):
                yield (a, b) if a <= b else (b, a)


@dataclass
class ScoreResult:
    tp: int = 0
    fn: int = 0
    fp: int = 0
    tn: int = 0
    must_link_total: int = 0
    cannot_link_total: int = 0
    must_link_uncovered: int = 0
    cannot_link_uncovered: int = 0
    # per-kind breakdown of the same counters
    by_kind: dict[str, dict[str, int]] = field(default_factory=dict)

    @property
    def precision(self) -> float:
        denom = self.tp + self.fp
        return self.tp / denom if denom else 1.0

    @property
    def recall(self) -> float:
        denom = self.tp + self.fn
        return self.tp / denom if denom else 1.0

    @property
    def f1(self) -> float:
        p, r = self.precision, self.recall
        return 2 * p * r / (p + r) if (p + r) else 0.0

    @property
    def must_link_covered(self) -> int:
        return self.tp + self.fn

    @property
    def cannot_link_covered(self) -> int:
        return self.fp + self.tn

    def as_dict(self) -> dict:
        return {
            "precision": round(self.precision, 4),
            "recall": round(self.recall, 4),
            "f1": round(self.f1, 4),
            "tp": self.tp,
            "fn": self.fn,
            "fp": self.fp,
            "tn": self.tn,
            "must_link": {
                "total": self.must_link_total,
                "covered": self.must_link_covered,
                "uncovered": self.must_link_uncovered,
            },
            "cannot_link": {
                "total": self.cannot_link_total,
                "covered": self.cannot_link_covered,
                "uncovered": self.cannot_link_uncovered,
            },
            "by_kind": self.by_kind,
        }


def _kind_bucket(kind: str) -> str:
    """Coarse label for the by-kind breakdown, robust to detailed wording."""
    k = kind.lower()
    if "numbered" in k:
        return "numbered_family"
    if "conflation" in k or "cross-entity" in k:
        return "conflation"
    return kind or "other"


def load_gold(path: str | Path) -> list[GoldCase]:
    """Load ground-truth cases from the partition JSONL format.

    Recognizes two record shapes:

    * ``clusters``: ``{label: [curie, ...]}`` — a cross-entity conflation case.
    * ``members``: ``{key: curie, ...}`` — a numbered family; every member is a
      distinct entity, so each becomes a singleton group.
    """
    cases: list[GoldCase] = []
    with open(path) as fh:
        for line_num, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            case = rec.get("case", f"case-{line_num}")
            kind = rec.get("kind", "")
            groups: list[frozenset[str]] = []
            if "clusters" in rec:
                for members in rec["clusters"].values():
                    if members:
                        groups.append(frozenset(members))
            elif "members" in rec:
                for curie in rec["members"].values():
                    groups.append(frozenset([curie]))
            else:
                raise ValueError(f"{path}:{line_num}: record has neither 'clusters' nor 'members': {case}")
            cases.append(
                GoldCase(
                    case=case,
                    kind=kind,
                    groups=tuple(groups),
                    needs_review=frozenset(rec.get("needs_review", [])),
                )
            )
    return cases


def membership_from_clusters(clusters: Iterable[Iterable[str]]) -> dict[str, int]:
    """Build a ``curie -> cluster_id`` map from an iterable of clusters."""
    membership: dict[str, int] = {}
    for cid, cluster in enumerate(clusters):
        for curie in cluster:
            membership[curie] = cid
    return membership


def load_membership(path: str | Path) -> dict[str, object]:
    """Load a membership map from JSON (``{curie: cluster}``) or JSONL rows
    (``{"curie": ..., "cluster": ...}``)."""
    path = Path(path)
    text = path.read_text().strip()
    if not text:
        return {}
    if text[0] == "{":
        return json.loads(text)
    membership: dict[str, object] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        row = json.loads(line)
        membership[row["curie"]] = row["cluster"]
    return membership


def score(gold: Iterable[GoldCase], membership: Membership) -> ScoreResult:
    """Score a membership map against ground-truth cases.

    Deduplicates pairs across cases; if a pair is asserted both must_link and
    cannot_link (a ground-truth contradiction), must_link wins and a warning is
    logged.
    """
    must: dict[tuple[str, str], str] = {}  # pair -> kind bucket
    cannot: dict[tuple[str, str], str] = {}

    for gc in gold:
        bucket = _kind_bucket(gc.kind)
        for pair in gc.must_link_pairs():
            must.setdefault(pair, bucket)
        for pair in gc.cannot_link_pairs():
            cannot.setdefault(pair, bucket)

    contradictions = must.keys() & cannot.keys()
    if contradictions:
        logging.warning(
            "Ground truth contradiction: %d pair(s) are both must_link and cannot_link; "
            "treating as must_link. Examples: %s",
            len(contradictions),
            sorted(contradictions)[:5],
        )
        for pair in contradictions:
            cannot.pop(pair, None)

    res = ScoreResult()
    res.must_link_total = len(must)
    res.cannot_link_total = len(cannot)

    def bump(bucket: str, counter: str) -> None:
        kd = res.by_kind.setdefault(bucket, {"tp": 0, "fn": 0, "fp": 0, "tn": 0})
        kd[counter] += 1

    def same_cluster(a: str, b: str) -> bool | None:
        if a not in membership or b not in membership:
            return None
        return membership[a] == membership[b]

    for (a, b), bucket in must.items():
        sc = same_cluster(a, b)
        if sc is None:
            res.must_link_uncovered += 1
        elif sc:
            res.tp += 1
            bump(bucket, "tp")
        else:
            res.fn += 1
            bump(bucket, "fn")

    for (a, b), bucket in cannot.items():
        sc = same_cluster(a, b)
        if sc is None:
            res.cannot_link_uncovered += 1
        elif sc:
            res.fp += 1
            bump(bucket, "fp")
        else:
            res.tn += 1
            bump(bucket, "tn")

    return res


def score_files(gold_path: str | Path, membership: str | Path | Membership) -> ScoreResult:
    gold = load_gold(gold_path)
    mem = membership if isinstance(membership, Mapping) else load_membership(membership)
    return score(gold, mem)


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Score an ER clustering against ground truth.")
    parser.add_argument("--gold", required=True, help="Path to ground-truth partition JSONL.")
    parser.add_argument(
        "--membership",
        required=True,
        help="Path to membership map (JSON {curie: cluster} or JSONL {curie, cluster}).",
    )
    parser.add_argument("--json", action="store_true", help="Emit the full result as JSON.")
    args = parser.parse_args(argv)

    result = score_files(args.gold, args.membership)
    if args.json:
        print(json.dumps(result.as_dict(), indent=2))
    else:
        print(
            f"precision={result.precision:.4f}  recall={result.recall:.4f}  f1={result.f1:.4f}\n"
            f"  must_link:   {result.must_link_covered} covered "
            f"({result.tp} same / {result.fn} apart), {result.must_link_uncovered} uncovered\n"
            f"  cannot_link: {result.cannot_link_covered} covered "
            f"({result.tn} apart / {result.fp} same), {result.cannot_link_uncovered} uncovered"
        )
        for bucket, counts in sorted(result.by_kind.items()):
            tp, fn, fp, tn = counts["tp"], counts["fn"], counts["fp"], counts["tn"]
            prec = tp / (tp + fp) if (tp + fp) else 1.0
            rec = tp / (tp + fn) if (tp + fn) else 1.0
            print(f"  [{bucket}] precision={prec:.4f} recall={rec:.4f} (tp={tp} fn={fn} fp={fp} tn={tn})")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
