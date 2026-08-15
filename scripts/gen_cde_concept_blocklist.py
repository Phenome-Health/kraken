"""Regenerate src/kraken/harmonizers/cde_concept_blocklist.py.

The CDE harmonizer emits `assesses`/`related_to` edges to NCIt/UMLS concepts drawn from the CDE
concept fields (objectClass/property/dataElementConcept) and from the constituents of compound
permissible-value codes. Many of those concepts are generic/grammatical glue (numbers, prepositions,
comparators, degree & temporal words, response-glue like Unknown/Other, phrase fragments) that make
poor edge targets. This script surfaces them automatically and a human curates the exceptions.

Method: collect every candidate concept CURIE, look each up in the integrated kraken node set, and
treat those typed ONLY as biolink:NamedThing as junk candidates (kraken couldn't type them). A human
keeps the substantive ones via RESCUE_IDS (concept-field rescues, by id) and RESCUE_ANSWER_NAMES
(answer-constituent rescues, by kraken name). Everything else is blocked. All id->name pairings in
the output come verbatim from the kraken node file -- none are hand-authored.

Re-run whenever the CDE export or the kraken node set changes. Requires ripgrep (falls back to grep).

    python scripts/gen_cde_concept_blocklist.py
"""

import json
import re
import shutil
import subprocess
from pathlib import Path

# Inputs: the raw CDE export and a built integrated kraken node file (for concept categories/names).
CDE_EXPORT = Path("/Volumes/AmySSD/kraken-data/input_data/cdes/SearchExport.json")
KRAKEN_NODES = Path("/Volumes/AmySSD/kraken-data/artifacts/integrated/kraken_nodes_2.0.1.jsonl")
OUT = Path(__file__).resolve().parents[1] / "src/kraken/harmonizers/cde_concept_blocklist.py"

NCIT_UMLS = {"nci thesaurus": "NCIT", "ncit": "NCIT", "nci": "NCIT", "umls": "UMLS"}

# Concept-field rescues (by id): substantive concepts under-typed as NamedThing in a SPOKE-less build.
RESCUE_IDS = {
    "NCIT:C25729",
    "NCIT:C71739",
    "NCIT:C189175",
    "NCIT:C25172",
    "NCIT:C17357",
    "NCIT:C179559",
    "NCIT:C164027",
    "NCIT:C164634",
    "NCIT:C21055",
    "NCIT:C16735",
    "NCIT:C38276",
    "NCIT:C38114",
    "NCIT:C28421",
    "NCIT:C124436",
    "NCIT:C83044",
    "NCIT:C83431",
    "NCIT:C178503",
    "NCIT:C178504",
    "NCIT:C124444",
    "NCIT:C25208",
}
# Answer-constituent rescues (by kraken node name, case-insensitive): meaningful response concepts.
RESCUE_ANSWER_NAMES = {
    s.lower()
    for s in {
        "Body Temperature",
        "Symptomatic",
        "Presymptomatic",
        "Transportation",
        "public transportation",
        "Sick Leave",
        "Paid Hours",
        "Unpaid Hours",
        "Hourly Employment",
        "Leave from Employment",
        "residence",
        "Home",
        "health status",
        "informed consent",
        "General Equivalency Diploma Completion",
        "Assent",
        # Answer values that also appear as compound fragments -- keep them so their standalone-answer
        # edges aren't dropped (frequency/severity/rating scales), consistent with un-blocked siblings.
        "monthly",
        "mild",
        "Moderate",
        "good",
    }
}


def collect_candidate_curies() -> set[str]:
    data = json.load(open(CDE_EXPORT))
    candidates: set[str] = set()
    for r in data:
        for field in ("objectClass", "property", "dataElementConcept"):
            for c in r.get(field, {}).get("concepts", []):
                origin, oid = (c.get("origin") or "").strip().lower(), c.get("originId")
                if origin in NCIT_UMLS and oid and ":" not in str(oid) and " " not in str(oid):
                    candidates.add(f"{NCIT_UMLS[origin]}:{oid}")
        for pv in r.get("valueDomain", {}).get("permissibleValues", []):
            for slot in ("conceptId", "valueMeaningCode"):
                v = pv.get(slot)
                if v and ":" in str(v):  # compound answer code -> split into constituents
                    for part in re.split(r"[:; ]+", str(v)):
                        m = re.match(r"^(C\d+)$", part.strip())
                        if m:
                            candidates.add("NCIT:" + m.group(1))
    return candidates


def lookup_categories_and_names(candidates: set[str]) -> tuple[dict, dict]:
    """Grep the (large) node file once for the candidate CURIEs; parse only matching lines."""
    cand_file = Path("/tmp/cde_blk_candidates.txt")
    cand_file.write_text("\n".join(sorted(candidates)))
    tool = shutil.which("rg") or shutil.which("grep")
    matches = subprocess.run(
        [tool, "-F", "-f", str(cand_file), str(KRAKEN_NODES)], capture_output=True, text=True
    ).stdout
    cat: dict[str, set[str]] = {}
    name: dict[str, str] = {}
    for line in matches.splitlines():
        try:
            n = json.loads(line)
        except json.JSONDecodeError:
            continue
        for t in (set(n.get("equivalent_ids", [])) | {n.get("id")}) & candidates:
            cat.setdefault(t, set()).update(n.get("categories", []))
            if n.get("name"):
                name.setdefault(t, n["name"])
    return cat, name


def main() -> None:
    candidates = collect_candidate_curies()
    cat, name = lookup_categories_and_names(candidates)
    namedthing_only = {t for t in candidates if cat.get(t) == {"biolink:NamedThing"}}
    rescued = {t for t in namedthing_only if t in RESCUE_IDS or name.get(t, "").lower() in RESCUE_ANSWER_NAMES}
    block = sorted(namedthing_only - rescued, key=lambda t: (name.get(t) or "").lower())

    header = '''"""
Curated blocklist of NCIt/UMLS concepts too generic/grammatical to be meaningful edge targets
(numbers, prepositions, comparators, degree & temporal glue, response-glue like Unknown/Other,
and phrase fragments). Covers both the concept fields (objectClass/property/dataElementConcept ->
assesses edges) and the constituents of compound permissible-value codes (-> related_to edges).

Generated by scripts/gen_cde_concept_blocklist.py -- do not edit by hand. Candidates are the concepts
typed ONLY as biolink:NamedThing in the kraken node set; substantive ones are kept via that script's
RESCUE_IDS / RESCUE_ANSWER_NAMES. All names below are verbatim from the kraken node file.

The harmonizer skips both the edge AND the stub node for any concept in this set.
"""

CONCEPT_BLOCKLIST: frozenset[str] = frozenset({'''
    lines = [header]
    for t in block:
        entry = repr(t) + ","
        nm = (name.get(t) or "").strip()
        if len(nm) > 80:
            nm = nm[:77] + "..."
        lines.append(f"    {entry:17} # {nm}")
    lines.append("})")
    lines.append("")
    OUT.write_text("\n".join(lines))
    print(
        f"candidates={len(candidates)} namedthing_only={len(namedthing_only)} "
        f"rescued={len(rescued)} blocked={len(block)}"
    )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
