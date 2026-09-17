"""Per-identifier facts from Babel -- each id's own label, categories and taxa -- for entity resolution.

Babel is ingested as a source (harmonizers/babel.py) with one node per identifier, so its facts ARE its harmonized
nodes. Entity resolution records them here in stage 1 and consults them in later stages: they are the source of
truth for an id's category and taxon, ahead of anything a source or a prefix guess says (see build.py).

Tens of millions of ids, so the store is sqlite on disk rather than a dict. It is rebuilt from the build's Babel
nodes on every run and deleted afterwards, so it runs without a journal or fsyncs.
"""

from __future__ import annotations

import json
import sqlite3
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class IdFacts:
    """One identifier's own facts. ``categories`` and ``taxa`` may be empty."""

    label: str | None
    categories: tuple[str, ...]
    taxa: tuple[str, ...] = ()


class IdFactsStore:
    def __init__(self, path: str | Path):
        self._db = sqlite3.connect(str(path))
        self._db.execute("PRAGMA journal_mode = OFF")
        self._db.execute("PRAGMA synchronous = OFF")
        self._db.execute(
            "CREATE TABLE IF NOT EXISTS facts (curie TEXT PRIMARY KEY, label TEXT, categories TEXT, taxa TEXT)"
        )

    def record(self, facts: Iterable[tuple[str, IdFacts]]) -> None:
        self._db.executemany(
            "INSERT OR REPLACE INTO facts (curie, label, categories, taxa) VALUES (?, ?, ?, ?)",
            (
                (curie, info.label, json.dumps(list(info.categories)), json.dumps(list(info.taxa)))
                for curie, info in facts
            ),
        )
        self._db.commit()

    def knows(self, curie: str) -> bool:
        """Whether Babel has this identifier at all (cheaper than ``get`` -- no row is read)."""
        return self._db.execute("SELECT 1 FROM facts WHERE curie = ?", (curie,)).fetchone() is not None

    def get(self, curie: str) -> IdFacts | None:
        row = self._db.execute("SELECT label, categories, taxa FROM facts WHERE curie = ?", (curie,)).fetchone()
        if row is None:
            return None
        label, categories, taxa = row
        return IdFacts(label=label, categories=tuple(json.loads(categories)), taxa=tuple(json.loads(taxa)))

    def resolve(self, curies: Iterable[str]) -> dict[str, IdFacts]:
        """Facts for each of ``curies`` that has any; ids Babel doesn't know are simply absent."""
        result: dict[str, IdFacts] = {}
        for curie in dict.fromkeys(curies):
            info = self.get(curie)
            if info is not None:
                result[curie] = info
        return result

    def close(self) -> None:
        self._db.commit()
        self._db.close()
