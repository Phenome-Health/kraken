"""Entity resolution's debug artifact: what every id is, where it went, and the evidence that put it there.

ER works through temp files -- Babel's cliques, every piece of match evidence, each id's names and types -- that
used to be deleted when it finished. This keeps the parts needed to answer "why is this id in this node?" after
the fact, in one SQLite file next to the build (``<integrated debug dir>/er_debug_<version>.sqlite``):

* ``ids``: one row per id in the build: the node it ended up in, Babel's own name / categories / taxon for it, the
  sources that provided it, and the name each source gave it.
* ``evidence``: every pair of ids with match evidence OTHER than a shared Babel clique -- a source's equivalence
  list, a same_as / close_match edge, a name match, an alias -- with its accumulated weight, whether it reached tau,
  and each kind of evidence behind it. (Babel pairs are left out because ``babel_cliques`` holds them far more
  compactly: two ids share a Babel clique iff they share a hub.)
* ``babel_cliques``: each Babel clique member's hub and clique size.

Read it with ``python -m kraken.entity_resolution.inspect_node CURIE``.
"""

from __future__ import annotations

import json
import sqlite3
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

BATCH = 200_000
DEBUG_DB_TEMPLATE = "er_debug_{version}.sqlite"

SCHEMA = """
CREATE TABLE ids (
    id TEXT PRIMARY KEY, node TEXT NOT NULL, babel_name TEXT, babel_categories TEXT, taxon TEXT,
    categories TEXT, provided_by TEXT, source_names TEXT
);
CREATE TABLE evidence (a TEXT NOT NULL, b TEXT NOT NULL, weight REAL NOT NULL, above_tau INTEGER NOT NULL,
    kinds TEXT NOT NULL);
CREATE TABLE babel_cliques (id TEXT PRIMARY KEY, hub TEXT NOT NULL, size INTEGER NOT NULL);
"""
INDEXES = """
CREATE INDEX ids_node ON ids(node);
CREATE INDEX evidence_a ON evidence(a);
CREATE INDEX evidence_b ON evidence(b);
CREATE INDEX babel_cliques_hub ON babel_cliques(hub);
"""


def debug_db_path(debug_dir: Path, version: str | None) -> Path:
    return Path(debug_dir) / DEBUG_DB_TEMPLATE.format(version=version or "build")


class DebugDbWriter:
    """Batched writer, used by the build. Indexes are built once, at ``close``."""

    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.unlink(missing_ok=True)
        self.path = path
        self._db = sqlite3.connect(path)
        self._db.executescript("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;" + SCHEMA)
        self._pending: dict[str, list[tuple]] = {"ids": [], "evidence": [], "babel_cliques": []}

    def _add(self, table: str, row: tuple) -> None:
        rows = self._pending[table]
        rows.append(row)
        if len(rows) >= BATCH:
            self._flush(table)

    def _flush(self, table: str) -> None:
        rows = self._pending[table]
        if rows:
            marks = ",".join("?" * len(rows[0]))
            self._db.executemany(f"INSERT OR REPLACE INTO {table} VALUES ({marks})", rows)
            self._db.commit()
            rows.clear()

    def add_id(
        self,
        curie: str,
        node: str,
        *,
        babel_name: str | None,
        babel_categories: Iterable[str],
        taxon: str | None,
        categories: Iterable[str],
        provided_by: Iterable[str],
        source_names: list[tuple[str, str]],
    ) -> None:
        self._add(
            "ids",
            (
                curie,
                node,
                babel_name,
                ",".join(babel_categories),
                taxon,
                ",".join(categories),
                ",".join(sorted(provided_by)),
                json.dumps(source_names) if source_names else None,
            ),
        )

    def add_evidence(self, a: str, b: str, weight: float, above_tau: bool, kinds: dict[str, float]) -> None:
        summary = ",".join(f"{kind}={weight:g}" for kind, weight in sorted(kinds.items()))
        self._add("evidence", (a, b, round(weight, 4), int(above_tau), summary))

    def add_babel_clique_member(self, curie: str, hub: str, size: int) -> None:
        self._add("babel_cliques", (curie, hub, size))

    def close(self) -> None:
        for table in self._pending:
            self._flush(table)
        self._db.executescript(INDEXES)
        self._db.commit()
        self._db.close()


@dataclass
class IdInfo:
    id: str
    node: str
    babel_name: str | None
    babel_categories: list[str]
    taxon: str | None
    categories: list[str]
    provided_by: list[str]
    source_names: list[tuple[str, str]]
    babel_hub: str | None = None
    babel_clique_size: int = 0


@dataclass
class Evidence:
    a: str
    b: str
    weight: float
    above_tau: bool
    kinds: str


class DebugDb:
    """Read side, for inspecting a build."""

    def __init__(self, path: Path):
        if not Path(path).exists():
            raise FileNotFoundError(f"no ER debug database at {path} (it is written by builds from this version on)")
        self._db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)

    def node_of(self, curie: str) -> str | None:
        row = self._db.execute("SELECT node FROM ids WHERE id = ?", (curie,)).fetchone()
        return row[0] if row else None

    def members(self, node: str) -> list[IdInfo]:
        rows = self._db.execute(
            "SELECT i.*, c.hub, c.size FROM ids i LEFT JOIN babel_cliques c ON c.id = i.id WHERE i.node = ? "
            "ORDER BY i.id",
            (node,),
        ).fetchall()
        return [
            IdInfo(
                id=r[0],
                node=r[1],
                babel_name=r[2],
                babel_categories=[c for c in (r[3] or "").split(",") if c],
                taxon=r[4],
                categories=[c for c in (r[5] or "").split(",") if c],
                provided_by=[p for p in (r[6] or "").split(",") if p],
                source_names=[tuple(pair) for pair in json.loads(r[7])] if r[7] else [],
                babel_hub=r[8],
                babel_clique_size=r[9] or 0,
            )
            for r in rows
        ]

    def evidence_touching(self, curies: Iterable[str]) -> Iterator[Evidence]:
        curies = sorted(set(curies))
        for start in range(0, len(curies), 500):
            chunk = curies[start : start + 500]
            marks = ",".join("?" * len(chunk))
            query = f"SELECT * FROM evidence WHERE a IN ({marks}) UNION SELECT * FROM evidence WHERE b IN ({marks})"
            for row in self._db.execute(query, chunk + chunk):
                yield Evidence(row[0], row[1], row[2], bool(row[3]), row[4])

    def babel_clique(self, hub: str) -> list[str]:
        return [r[0] for r in self._db.execute("SELECT id FROM babel_cliques WHERE hub = ? ORDER BY id", (hub,))]
