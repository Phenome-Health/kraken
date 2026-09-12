"""Every curie KRAKEN ingests goes through biomapper2 (``BaseHarmonizer.normalize_curie``).

That matters most for the ids an aggregator records as an edge's ORIGINAL endpoints: entity
resolution treats those as real ids -- they seed the match graph and edges are remapped through them
-- so they must spell their vocabulary the way the rest of KRAKEN does. Translator writes
``Ensembl:ENSG00000099864`` while the node set uses ``ENSEMBL:``; left alone those are two different
ids and the edge finds no node.

biomapper2's own lookups are tested in biomapper2 (and exercised here against real source data);
these cover the decisions this layer makes around it -- above all that an id is never dropped.
"""

from collections import defaultdict

import pytest

from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import ORIGINAL_OBJECT_ATTR, ORIGINAL_SUBJECT_ATTR


class _Harmonizer(BaseHarmonizer):
    source_infores = "infores:test"


class _StubNormalizer:
    """Resolves only what it is told about; everything else is an unrecognized vocab.

    ``curies`` maps a (vocab, local_id) pair to the curie biomapper2 would return. A vocab listed in
    ``known_vocabs`` but with no mapping for the id is treated as a recognized vocab with an invalid
    local id -- the other way normalization can fail, which is reported separately.
    """

    def __init__(self, curies: dict[tuple[str, str], str], known_vocabs: set[str] | None = None):
        self.curies = curies
        self.known_vocabs = known_vocabs or {vocab for vocab, _ in curies}

    def get_curies(self, local_ids_dict, **_kwargs):
        ((vocab, local_id),) = local_ids_dict.items()
        prefix = vocab if isinstance(vocab, str) else vocab[0]
        resolved = self.curies.get((prefix, local_id))
        if resolved:
            return {resolved: ""}, {}, set()
        if prefix in self.known_vocabs:
            return {}, {prefix: [local_id]}, set()  # known vocab, bad id
        return {}, {}, {prefix}  # unrecognized vocab


def _harmonizer(curies=None, known_vocabs=None) -> BaseHarmonizer:
    """A bare harmonizer (BaseHarmonizer.__init__ builds a Biolink toolkit and a real Normalizer)."""
    instance = object.__new__(_Harmonizer)
    instance.normalizer = _StubNormalizer(curies or {}, known_vocabs)
    instance.normalized_id_map = {}
    instance.unrecognized_vocabs = set()
    instance.prefixes_with_invalid_ids = defaultdict(int)
    instance.unrecognized_vocab_prefixes = {}
    instance.invalid_id_prefixes = {}
    instance.invalid_curies = set()
    return instance


def test_prefix_is_standardized():
    h = _harmonizer({("Ensembl", "ENSG00000099864"): "ENSEMBL:ENSG00000099864"})
    assert h.normalize_curie("Ensembl:ENSG00000099864") == "ENSEMBL:ENSG00000099864"


def test_local_id_containing_colons_is_split_only_at_the_prefix():
    """An HGVS expression's local id contains colons; only the first one ends the prefix."""
    h = _harmonizer({("hgvs", "NC_000021.9:g.25840043C>G"): "HGVS:NC_000021.9:g.25840043C>G"})
    assert h.normalize_curie("hgvs:NC_000021.9:g.25840043C>G") == "HGVS:NC_000021.9:g.25840043C>G"


def test_unrecognized_vocabulary_is_kept_not_dropped():
    """The id survives untouched, and is tallied as a vocab biomapper2 should learn."""
    h = _harmonizer()
    assert h.normalize_curie("WEIRD:12345") == "WEIRD:12345"
    assert h.unrecognized_vocab_prefixes["WEIRD"]["count"] == 1
    assert h.unrecognized_vocab_prefixes["WEIRD"]["examples"] == ["WEIRD:12345"]
    assert not h.invalid_id_prefixes  # a different failure, reported separately


def test_invalid_local_id_is_kept_and_tallied_separately():
    """Vocab recognized, id rejected -- a different fix from 'add this vocab', so counted apart."""
    h = _harmonizer(known_vocabs={"DRUGBANK"})
    assert h.normalize_curie("DRUGBANK:not-an-id") == "DRUGBANK:not-an-id"
    assert h.invalid_id_prefixes["DRUGBANK"]["count"] == 1
    assert not h.unrecognized_vocab_prefixes


def test_tallies_count_distinct_curies_and_cap_examples():
    h = _harmonizer()
    for i in range(10):
        h.normalize_curie(f"WEIRD:{i}")
    h.normalize_curie("WEIRD:0")  # repeat -- already cached, must not double-count
    entry = h.unrecognized_vocab_prefixes["WEIRD"]
    assert entry["count"] == 10
    assert len(entry["examples"]) == 3


def test_non_curie_passes_through():
    h = _harmonizer()
    assert h.normalize_curie("not-a-curie") == "not-a-curie"
    assert not h.unrecognized_vocab_prefixes  # nothing to report; it isn't a curie at all


def test_result_is_cached():
    h = _harmonizer({("Ensembl", "E1"): "ENSEMBL:E1"})
    calls = []
    inner = h.normalizer.get_curies

    def counting(local_ids_dict, **kwargs):
        calls.append(local_ids_dict)
        return inner(local_ids_dict, **kwargs)

    h.normalizer.get_curies = counting
    assert h.normalize_curie("Ensembl:E1") == "ENSEMBL:E1"
    assert h.normalize_curie("Ensembl:E1") == "ENSEMBL:E1"
    assert len(calls) == 1  # second call served from the cache


def test_kegg_is_offered_its_sub_vocabularies():
    """Sources use a bare KEGG prefix where they mean KEGG.COMPOUND; biomapper2 picks."""
    h = _harmonizer({("kegg", "C00031"): "KEGG.COMPOUND:C00031"})
    assert h.normalize_curie("KEGG:C00031") == "KEGG.COMPOUND:C00031"


@pytest.mark.parametrize("attr", [ORIGINAL_SUBJECT_ATTR, ORIGINAL_OBJECT_ATTR])
def test_original_endpoints_are_normalized_in_place(attr):
    h = _harmonizer({("Ensembl", "E1"): "ENSEMBL:E1"})
    attributes = {attr: "Ensembl:E1", "p_value": 0.01}
    h._normalize_original_endpoints(attributes)
    assert attributes[attr] == "ENSEMBL:E1"
    assert attributes["p_value"] == 0.01  # everything else left alone


def test_missing_or_non_string_originals_are_ignored():
    h = _harmonizer()
    attributes = {ORIGINAL_SUBJECT_ATTR: None, "other": "x"}
    h._normalize_original_endpoints(attributes)
    assert attributes == {ORIGINAL_SUBJECT_ATTR: None, "other": "x"}
