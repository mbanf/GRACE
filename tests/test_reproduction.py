"""Reproduction of the published GRACE v1.0 results from the tutorial datasets.

Needs the unpacked datasets; set GRACE_DATA_DIR to the folder holding athaliana_datasets/ and
drosophila_datasets/. Skipped otherwise.
"""

import os
from pathlib import Path

import pytest

from grace.pipeline import ensemble_network, prepare
from grace.reference import published_models, published_network

DATA = Path(os.environ.get("GRACE_DATA_DIR", "data"))
SPECIES = [s for s in ("athaliana", "drosophila") if (DATA / f"{s}_datasets").exists()]
pytestmark = pytest.mark.skipif(not SPECIES, reason="set GRACE_DATA_DIR to the tutorial datasets")

# sizes of the inputs as computed by the R implementation (verified against R)
EXPECTED = {
    "athaliana": dict(links=10098, modules=560, tree_edges=1995, coreg_entries=7312, gold_standard=44,
                      cofunctional_evidence=10626, network=792,
                      validation={"AtRegNet": 1577, "SUBA_colocalization": 11452468}),
    "drosophila": dict(links=17772, modules=706, tree_edges=7457, coreg_entries=34066, gold_standard=10,
                       cofunctional_evidence=38423, network=7164,
                       validation={"HiC": 86632, "ChIP": 75499, "ImaGO_colocalization": 4504778}),
}


@pytest.fixture(scope="module", params=SPECIES)
def prepared(request):
    return request.param, prepare(request.param, DATA)


def test_inputs_match_r(prepared):
    species, p = prepared
    e = EXPECTED[species]
    assert len(p.network.links) == e["links"]
    assert p.network.modules.module.nunique() == e["modules"]
    assert len(p.network.edges) == e["tree_edges"]
    assert p.network.adj.nnz == e["coreg_entries"]
    assert len(p.gs_pairs) == e["gold_standard"]
    assert p.training.bp.nnz == e["cofunctional_evidence"]
    for k, n in e["validation"].items():
        assert len(p.dataset.validation[k]) == n


def test_published_models_reproduce_published_network(prepared):
    species, p = prepared
    net = ensemble_network(p, published_models(species))
    assert set(zip(net.TF, net.TG)) == published_network(species)
    assert len(net) == EXPECTED[species]["network"]
