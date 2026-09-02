"""Tests for the SRI Node Normalizer client (cache, inference, parsing)."""

from kraken.entity_resolution.sri_nodenorm import NodeNormClient, NormInfo, infer_category


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


def test_infer_category():
    assert infer_category("HGNC:2707") == "biolink:Gene"
    assert infer_category("UniProtKB:P12821") == "biolink:Protein"
    assert infer_category("CHEBI:1234") == "biolink:ChemicalEntity"
    assert infer_category("WEIRD:1") is None


def test_inference_is_backup_when_nn_unrecognized(tmp_path):
    # NN returns nothing -> prefix inference supplies the category as a BACKUP,
    # but the API was still consulted first (inference never pre-empts it).
    session = _FakeSession({})
    client = NodeNormClient(tmp_path / "c.sqlite", session=session)
    out = client.resolve(["HGNC:2707", "NCBITaxon:9606"])
    assert out["HGNC:2707"].categories == ("biolink:Gene",)
    assert out["NCBITaxon:9606"].categories == ("biolink:OrganismTaxon",)
    assert session.calls == 1  # one batch queried; inference only filled the gap
    client.close()


def test_nn_category_overrides_inference(tmp_path):
    # Even for an inferable prefix (HGNC->Gene), the normalizer's answer wins.
    payload = {
        "HGNC:2707": {
            "equivalent_identifiers": [{"identifier": "HGNC:2707", "label": "ACE", "type": ["biolink:Protein"]}]
        }
    }
    session = _FakeSession(payload)
    client = NodeNormClient(tmp_path / "c.sqlite", session=session)
    out = client.resolve(["HGNC:2707"])
    assert out["HGNC:2707"] == NormInfo(label="ACE", categories=("biolink:Protein",))
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
    parsed = NodeNormClient._parse_response(["FOO:1"], data)
    # per-member label (NOT the clique's id.label), and per-member type
    assert parsed["FOO:1"] == NormInfo(label="member label", categories=("biolink:Protein",))


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
