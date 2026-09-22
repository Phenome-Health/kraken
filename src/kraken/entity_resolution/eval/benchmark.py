"""Score a KRAKEN build's entity resolution against the benchmark cases.

Each case is a real concept, curated from the real data, partitioned into the distinct entities it involves:

    {"case": "atrial-fibrillation", "kind": "disease",
     "clusters": {"atrial fibrillation": ["MONDO:0004981", "HP:0005110", ...],
                  "atrial flutter": ["MONDO:0004940", ...]},
     "known_issue": "why it fails today, if it does",
     "flags": [{"ids": ["MP:0003705"], "note": "best guess: a mouse phenotype, but it names the disease"}],
     "notes": "how the case was built"}

A case PASSES when every id is in the build, each cluster's ids are in one node, and no two clusters share a node.
Cases marked ``known_issue`` document what the build still gets wrong: a known issue that passes is FIXED, and any
other case that fails is a REGRESSION. ``flags`` mark ids whose placement is a best guess, for review; they are
scored like the rest.

The build is read from its nodes JSONL (``--nodes``), or from a running Kestrel API (``--api``), which answers
``/canonicalize`` without reading an 11 GB file.

    uv run python -m kraken.entity_resolution.eval.benchmark --nodes kraken_nodes_2.1.1.jsonl
    uv run python -m kraken.entity_resolution.eval.benchmark --api https://kestrel.krakenkg.com/api
    uv run python -m kraken.entity_resolution.eval.benchmark --api http://localhost:9990/api --json > report.json

The same cases also give pairwise precision/recall (see ``scorer``), over every must-link pair (two ids of one
cluster) and cannot-link pair (ids of two clusters of one case).
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

from kraken.entity_resolution.eval.scorer import DEFAULT_GROUND_TRUTH_PATH, load_gold, score
from kraken.utils.constants import PROJECT_ROOT

DEFAULT_BENCHMARK_PATH = PROJECT_ROOT / "config" / "entity_resolution" / "er_benchmark.jsonl"
API_BATCH = 500


@dataclass
class BenchmarkCase:
    case: str
    kind: str
    clusters: dict[str, list[str]]
    known_issue: str | None = None
    flags: list[dict] = field(default_factory=list)
    notes: str | None = None

    @property
    def ids(self) -> set[str]:
        return {curie for members in self.clusters.values() for curie in members}

    @property
    def flagged_ids(self) -> set[str]:
        return {curie for flag in self.flags for curie in flag.get("ids", [])}


@dataclass
class CaseResult:
    case: BenchmarkCase
    missing: list[str]
    # cluster label -> {node id: [ids]} for every cluster spread over more than one node
    splits: dict[str, dict[str, list[str]]]
    # node id -> {cluster label: [ids]} for every node holding ids of more than one cluster
    conflations: dict[str, dict[str, list[str]]]

    @property
    def passed(self) -> bool:
        return not (self.missing or self.splits or self.conflations)

    @property
    def outcome(self) -> str:
        if self.case.known_issue:
            return "fixed" if self.passed else "known issue"
        return "pass" if self.passed else "REGRESSION"

    def as_dict(self) -> dict:
        return {
            "case": self.case.case,
            "kind": self.case.kind,
            "outcome": self.outcome,
            "known_issue": self.case.known_issue,
            "missing": self.missing,
            "splits": self.splits,
            "conflations": self.conflations,
            "flagged_ids": sorted(self.case.flagged_ids),
        }


def load_cases(path: str | Path) -> list[BenchmarkCase]:
    cases = []
    with open(path) as fh:
        for line_num, line in enumerate(fh, 1):
            if not line.strip():
                continue
            rec = json.loads(line)
            if "clusters" not in rec:
                raise ValueError(f"{path}:{line_num}: a benchmark case needs 'clusters': {rec.get('case')}")
            cases.append(
                BenchmarkCase(
                    case=rec["case"],
                    kind=rec.get("kind", ""),
                    clusters=rec["clusters"],
                    known_issue=rec.get("known_issue"),
                    flags=rec.get("flags", []),
                    notes=rec.get("notes"),
                )
            )
    names = [c.case for c in cases]
    duplicates = {n for n in names if names.count(n) > 1}
    if duplicates:
        raise ValueError(f"{path}: duplicate case names {sorted(duplicates)}")
    return cases


def membership_from_nodes_file(path: str | Path, wanted: set[str]) -> dict[str, str]:
    """``curie -> node id`` for every wanted curie in a nodes JSONL (streamed; only the wanted ids are kept)."""
    membership: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            node = json.loads(line)
            for curie in node.get("equivalent_ids") or [node["id"]]:
                if curie in wanted:
                    membership[curie] = node["id"]
    return membership


def membership_from_api(base_url: str, wanted: set[str]) -> dict[str, str]:
    """``curie -> node id`` from a Kestrel API's ``/canonicalize``; ids it doesn't know are left out."""
    membership: dict[str, str] = {}
    curies = sorted(wanted)
    for start in range(0, len(curies), API_BATCH):
        body = json.dumps({"curies": curies[start : start + API_BATCH]}).encode()
        request = urllib.request.Request(
            f"{base_url.rstrip('/')}/canonicalize", data=body, headers={"content-type": "application/json"}
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            for curie, node_id in json.load(response).items():
                if node_id:
                    membership[curie] = node_id
    return membership


def evaluate(case: BenchmarkCase, membership: Mapping[str, str]) -> CaseResult:
    missing = sorted(curie for curie in case.ids if curie not in membership)
    nodes_of: dict[str, dict[str, list[str]]] = {}  # cluster label -> node -> ids
    labels_on: dict[str, dict[str, list[str]]] = defaultdict(dict)  # node -> cluster label -> ids
    for label, members in case.clusters.items():
        by_node: dict[str, list[str]] = defaultdict(list)
        for curie in members:
            if curie in membership:
                by_node[membership[curie]].append(curie)
        nodes_of[label] = dict(by_node)
        for node, ids in by_node.items():
            labels_on[node][label] = ids
    splits = {label: by_node for label, by_node in nodes_of.items() if len(by_node) > 1}
    conflations = {node: labels for node, labels in labels_on.items() if len(labels) > 1}
    return CaseResult(case, missing, splits, conflations)


def run(cases: Iterable[BenchmarkCase], membership: Mapping[str, str]) -> list[CaseResult]:
    return [evaluate(case, membership) for case in cases]


def _pairwise(case_files: Iterable[Path], membership: Mapping[str, str]) -> dict:
    gold = [g for path in case_files for g in load_gold(path)]
    return score(gold, membership).as_dict()


def _print_report(results: list[CaseResult], pairwise: dict, out=sys.stdout) -> None:
    counts = defaultdict(int)
    for result in results:
        counts[result.outcome] += 1
    by_kind: dict[str, list[CaseResult]] = defaultdict(list)
    for result in results:
        by_kind[result.case.kind or "other"].append(result)

    def line(text: str = "") -> None:
        print(text, file=out)

    line(f"ER benchmark: {len(results)} cases")
    line(
        f"  pass {counts['pass']}   REGRESSION {counts['REGRESSION']}   "
        f"known issue {counts['known issue']}   fixed {counts['fixed']}"
    )
    line(
        f"  pairwise precision={pairwise['precision']:.4f}  recall={pairwise['recall']:.4f}  f1={pairwise['f1']:.4f}"
        f"  (must-link {pairwise['tp']}/{pairwise['tp'] + pairwise['fn']} together,"
        f" cannot-link {pairwise['tn']}/{pairwise['tn'] + pairwise['fp']} apart)"
    )
    line()
    line("By kind:")
    for kind, kind_results in sorted(by_kind.items()):
        passed = sum(r.passed for r in kind_results)
        line(f"  {kind:<22} {passed}/{len(kind_results)} passing")
    for heading, outcome in (("REGRESSIONS", "REGRESSION"), ("FIXED", "fixed"), ("KNOWN ISSUES", "known issue")):
        chosen = [r for r in results if r.outcome == outcome]
        if not chosen:
            continue
        line()
        line(f"{heading} ({len(chosen)}):")
        for result in chosen:
            line(f"  {result.case.case}" + (f"  -- {result.case.known_issue}" if result.case.known_issue else ""))
            if result.missing:
                line(f"      missing from the build: {_ids(result.missing)}")
            for label, by_node in result.splits.items():
                # the node holding most of the cluster is where it "is"; the rest are the strays
                main = max(by_node, key=lambda node: len(by_node[node]))
                strays = "; ".join(f"{node}: {_ids(ids)}" for node, ids in sorted(by_node.items()) if node != main)
                line(f"      split '{label}' ({len(by_node[main])} ids in {main}) -- also in {strays}")
            for node, labels in result.conflations.items():
                parts = "; ".join(f"'{label}' ({_ids(ids)})" for label, ids in sorted(labels.items()))
                line(f"      merged in {node}: {parts}")
    flagged = [r for r in results if r.case.flagged_ids]
    if flagged:
        line()
        line(f"{len(flagged)} cases have flagged (best-guess) ids; list them with --flags")
    line()
    line("To see every id of a node and why it's there: python -m kraken.entity_resolution.inspect_node <any id>")


def _ids(ids: list[str], shown: int = 6) -> str:
    return ", ".join(ids[:shown]) + (f", ... ({len(ids)} ids)" if len(ids) > shown else "")


def _print_flags(cases: Iterable[BenchmarkCase], out=sys.stdout) -> None:
    for case in cases:
        for flag in case.flags:
            print(f"{case.case}: {_ids(flag.get('ids', []))}\n    {flag.get('note', '')}", file=out)


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Score a KRAKEN build against the entity-resolution benchmark.")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--nodes", help="the build's nodes JSONL")
    source.add_argument("--api", help="a Kestrel API base url, e.g. https://kestrel.krakenkg.com/api")
    parser.add_argument("--cases", action="append", help=f"case file(s) (default: {DEFAULT_BENCHMARK_PATH.name})")
    parser.add_argument(
        "--with-seed", action="store_true", help=f"also score the pairwise metrics on {DEFAULT_GROUND_TRUTH_PATH.name}"
    )
    parser.add_argument("--json", action="store_true", help="emit the full report as JSON")
    parser.add_argument("--flags", action="store_true", help="list every flagged (best-guess) id for review, and exit")
    args = parser.parse_args(argv)

    case_files = [Path(p) for p in (args.cases or [DEFAULT_BENCHMARK_PATH])]
    cases = [case for path in case_files for case in load_cases(path)]
    if args.flags:
        _print_flags(cases)
        return 0
    if not (args.nodes or args.api):
        parser.error("one of --nodes or --api is required")
    pairwise_files = case_files + ([DEFAULT_GROUND_TRUTH_PATH] if args.with_seed else [])
    wanted = {curie for case in cases for curie in case.ids}
    wanted |= {curie for g in load_gold(DEFAULT_GROUND_TRUTH_PATH) for group in g.groups for curie in group}
    if args.nodes:
        membership = membership_from_nodes_file(args.nodes, wanted)
    else:
        membership = membership_from_api(args.api, wanted)

    results = run(cases, membership)
    pairwise = _pairwise(pairwise_files, membership)
    if args.json:
        print(json.dumps({"pairwise": pairwise, "cases": [r.as_dict() for r in results]}, indent=1))
    else:
        _print_report(results, pairwise)
    return 1 if any(r.outcome == "REGRESSION" for r in results) else 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
