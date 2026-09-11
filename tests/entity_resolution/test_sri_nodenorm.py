"""Tests for the SRI Node Normalizer client (cache, inference, parsing, retries)."""

import requests

from kraken.entity_resolution.sri_nodenorm import NodeNormClient, NormInfo, infer_category, infer_taxon


class _FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        pass

    def json(self):
        return self._payload


class _FakeSession:
    def __init__(self, payload):
        self.payload = payload
        self.calls = 0

    def post(self, url, json, timeout):  # noqa: A002 - mirror requests signature
        self.calls += 1
        self.last_curies = json["curies"]
        return _FakeResponse(self.payload)


class _FlakySession:
    """Raises a connection error for the first ``fail_times`` calls, then succeeds."""

    def __init__(self, payload, fail_times):
        self.payload = payload
        self.fail_times = fail_times
        self.calls = 0

    def post(self, url, json, timeout):  # noqa: A002
        self.calls += 1
        if self.calls <= self.fail_times:
            raise requests.ConnectionError("boom")
        return _FakeResponse(self.payload)


def test_retries_then_succeeds(tmp_path):
    payload = {
        "FOO:1": {
            "equivalent_identifiers": [{"identifier": "FOO:1", "label": "Foo", "type": ["biolink:SmallMolecule"]}]
        }
    }
    session = _FlakySession(payload, fail_times=2)
    client = NodeNormClient(tmp_path / "c.sqlite", session=session, retry_backoff=0.0, max_retries=3)
    out = client.resolve(["FOO:1"])
    assert out["FOO:1"].categories == ("biolink:SmallMolecule",)
    assert session.calls == 3  # 2 failures + 1 success
    client.close()


def test_failed_batch_is_not_cached_and_retried_next_run(tmp_path):
    # A batch that fails after all retries must NOT be cached (else a transient blip
    # poisons the cache); it should be re-queried on the next run.
    session = _FlakySession({}, fail_times=99)  # always fails
    client = NodeNormClient(tmp_path / "c.sqlite", session=session, retry_backoff=0.0, max_retries=2)
    out = client.resolve(["FOO:9"])
    assert out["FOO:9"] == NormInfo(label=None, categories=())  # fallback for this run
    calls_first = session.calls  # 1 + 2 retries = 3
    assert calls_first == 3
    client.resolve(["FOO:9"])  # not cached -> re-queried
    assert session.calls > calls_first
    client.close()


def test_infer_category():
    assert infer_category("HGNC:2707") == "biolink:Gene"
    assert infer_category("UniProtKB:P12821") == "biolink:Protein"
    assert infer_category("CHEBI:1234") == "biolink:ChemicalEntity"
    assert infer_category("WEIRD:1") is None


def test_resolve_does_not_backfill_category_or_taxon_for_unrecognized(tmp_path):
    # A prefix guess (category OR taxon) must NEVER pre-empt a source-derived value, so
    # resolve() itself does NOT backfill either for an NN-unrecognized id -- it returns
    # empty facts. The CALLER (build.py) applies infer_category()/infer_taxon() only AFTER
    # source values.
    session = _FakeSession({})
    client = NodeNormClient(tmp_path / "c.sqlite", session=session)
    out = client.resolve(["HGNC:2707", "NCBITaxon:9606"])
    assert out["HGNC:2707"].categories == ()  # no category backup in resolve()
    assert out["HGNC:2707"].taxa == ()  # no taxon backup in resolve() either
    assert out["NCBITaxon:9606"].categories == ()
    assert session.calls == 1  # the API was still consulted first
    # The guesses remain available for the caller to use as a last resort:
    assert infer_category("HGNC:2707") == "biolink:Gene"
    assert infer_taxon("HGNC:2707") == "NCBITaxon:9606"
    client.close()


def test_nn_answer_is_returned_verbatim_without_prefix_backup(tmp_path):
    # The normalizer's own answer is returned as-is; resolve() adds no prefix backup, so a
    # taxon NN omitted is left empty here (the caller backfills it after source).
    payload = {
        "HGNC:2707": {
            "equivalent_identifiers": [{"identifier": "HGNC:2707", "label": "ACE", "type": ["biolink:Protein"]}]
        }
    }
    session = _FakeSession(payload)
    client = NodeNormClient(tmp_path / "c.sqlite", session=session)
    out = client.resolve(["HGNC:2707"])
    assert out["HGNC:2707"] == NormInfo(label="ACE", categories=("biolink:Protein",), taxa=())
    client.close()


def test_parse_response_uses_per_member_label_and_type():
    data = {
        "FOO:1": {
            "id": {"identifier": "FOO:1", "label": "clique preferred label"},
            "type": ["biolink:Gene"],
            "equivalent_identifiers": [
                {"identifier": "FOO:1", "label": "member label", "type": ["biolink:Protein"]},
                {"identifier": "FOO:2", "label": "other", "type": ["biolink:Gene"]},
            ],
        }
    }
    parsed = NodeNormClient._parse_response(data)
    # HARVEST: every clique member is returned (not just the queried id), each with
    # its per-member label/type (NOT the clique id.label) and the clique canonical.
    assert parsed["FOO:1"] == NormInfo(
        label="member label", categories=("biolink:Protein",), canonical="FOO:1"
    )
    assert parsed["FOO:2"] == NormInfo(label="other", categories=("biolink:Gene",), canonical="FOO:1")


def test_api_called_for_uninferable_and_cached(tmp_path):
    payload = {
        "FOO:1": {
            "equivalent_identifiers": [{"identifier": "FOO:1", "label": "Foo one", "type": ["biolink:SmallMolecule"]}]
        }
    }
    session = _FakeSession(payload)
    client = NodeNormClient(tmp_path / "c.sqlite", session=session)
    out = client.resolve(["FOO:1"])
    assert out["FOO:1"] == NormInfo(label="Foo one", categories=("biolink:SmallMolecule",))
    assert session.calls == 1
    # second call served from cache -> no new API call
    out2 = client.resolve(["FOO:1"])
    assert out2["FOO:1"].categories == ("biolink:SmallMolecule",)
    assert session.calls == 1
    client.close()


def test_negative_cached(tmp_path):
    session = _FakeSession({})  # empty payload -> unresolved
    client = NodeNormClient(tmp_path / "c.sqlite", session=session)
    out = client.resolve(["FOO:99"])
    assert out["FOO:99"] == NormInfo(label=None, categories=())
    assert session.calls == 1
    client.resolve(["FOO:99"])  # negative is cached
    assert session.calls == 1
    client.close()


def test_harvest_dedups_clique_mates_and_exposes_cliques(tmp_path):
    # Querying ANY member returns the whole clique; clique-mates get cached, so a
    # later queued member is a cache hit (not a second API call), and the clique is
    # exposed via iter_cliques as the equivalence signal.
    payload = {
        "NCBIGene:7157": {
            "id": {"identifier": "NCBIGene:7157"},
            "type": ["biolink:Gene"],
            "taxa": ["NCBITaxon:9606"],
            "equivalent_identifiers": [
                {"identifier": "NCBIGene:7157", "label": "TP53", "type": ["biolink:Gene"], "taxa": ["NCBITaxon:9606"]},
                {"identifier": "HGNC:11998", "label": "TP53", "type": ["biolink:Gene"], "taxa": ["NCBITaxon:9606"]},
            ],
        }
    }
    session = _FakeSession(payload)
    client = NodeNormClient(tmp_path / "c.sqlite", session=session, batch_size=1)
    out = client.resolve(["NCBIGene:7157", "HGNC:11998"])
    assert session.calls == 1  # HGNC came from the gene's clique, not a second query
    assert out["HGNC:11998"].canonical == "NCBIGene:7157"
    assert out["HGNC:11998"].taxa == ("NCBITaxon:9606",)
    cliques = list(client.iter_cliques())
    assert len(cliques) == 1
    canonical, members = cliques[0]
    assert canonical == "NCBIGene:7157"
    assert set(members) == {"NCBIGene:7157", "HGNC:11998"}
    client.close()


def test_taxon_prefix_backup_is_the_callers_job_not_resolves(tmp_path):
    # resolve() does NOT backfill taxon: an NN-unrecognized id comes back untaxoned, and
    # the single-species prefix guess is exposed via infer_taxon() for the caller to apply
    # AFTER the source taxon (see build.py mg_taxon / node-taxon finalization).
    session = _FakeSession({})  # nothing recognized
    client = NodeNormClient(tmp_path / "c.sqlite", session=session)
    out = client.resolve(["MGI:98834"])  # MGI == mouse by construction
    assert out["MGI:98834"].taxa == ()  # resolve() adds no taxon backup
    assert infer_taxon("MGI:98834") == "NCBITaxon:10090"  # the guess is available to the caller
    client.close()
