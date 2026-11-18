"""Data-free tests of the GRACE core: exact MRF inference, network construction, co-regulation."""

import itertools

import numpy as np
import pandas as pd
import pytest

from grace.mrf import build_forest, predict_links
from grace.network import build_network
from grace.validation import coregulated_pairs


def brute_force_marginals(v, edges, alpha, gamma, lam):
    """P(link) per node by enumerating all states of a small Potts tree with clamping."""
    n = len(v)
    clamped = v > alpha
    p1 = np.zeros(n)
    z = 0.0
    for states in itertools.product([0, 1], repeat=n):   # 1 = link
        s = np.array(states)
        if np.any(clamped & (s == 0)):
            continue
        w = np.prod(np.where(s == 1, np.exp(v - gamma), 1.0))
        for a, b, wt in edges:
            if s[a] == s[b]:
                w *= np.exp(lam * wt)
        z += w
        p1 += w * s
    return p1 / z


def random_tree(rng, n):
    return [(k, int(rng.integers(0, k)), float(rng.uniform(0.1, 0.95))) for k in range(1, n)]


@pytest.mark.parametrize("seed", range(25))
def test_tree_inference_matches_brute_force(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(2, 11))
    v = rng.uniform(0.0, 0.03, n)
    edges = random_tree(rng, n)
    gamma = float(rng.uniform(0.005, 0.025))
    lam = float(rng.uniform(0.01, 3.0))
    forest = build_forest(n, np.zeros(n, int), np.arange(n), np.arange(1, n + 1),
                          np.zeros(len(edges), int), np.array([a + 1 for a, _, _ in edges]),
                          np.array([b + 1 for _, b, _ in edges]), np.array([w for *_, w in edges]))
    import scipy.sparse as sp
    a = np.array([[x, y] for x, y, _ in edges])
    adj = sp.csr_matrix((np.ones(2 * len(a)), (np.r_[a[:, 0], a[:, 1]], np.r_[a[:, 1], a[:, 0]])), shape=(n, n))
    out = predict_links(forest, v, adj.indptr, adj.indices, gamma, gamma, lam)
    p = brute_force_marginals(v, edges, gamma, gamma, lam)
    expected = np.where((v > gamma) | (p > 0.5), 1.0, 0.0)
    clear = np.abs(p - 0.5) > 1e-9                     # skip exact ties
    assert np.array_equal(out[clear], expected[clear])


def test_singletons_use_threshold_only():
    from grace.mrf import Forest

    v = np.array([0.02, 0.005])
    empty = Forest(np.zeros(0, int), np.full(2, -1), np.zeros(2), np.full(2, -1), np.zeros(0, int), 0)
    out = predict_links(empty, v, np.zeros(3, dtype=np.int32), np.zeros(0, dtype=np.int32), 0.01, 0.01, 1.0)
    assert out.tolist() == [1.0, 0.0]


def test_build_network_links_modules_and_trees():
    grn = pd.DataFrame({"TF": ["A", "A", "A", "B", "B"], "TG": ["x", "y", "z", "x", "w"],
                        "value": [0.9, 0.5, 0.4, 0.3, 0.0]})
    cof = pd.DataFrame({"row": ["x", "y", "x", "z", "y", "z"], "col": ["y", "x", "z", "x", "z", "y"],
                        "value": [0.8, 0.8, 0.2, 0.2, 0.6, 0.6]})
    net = build_network(grn, cof)
    assert list(net.links.TF) == ["A", "A", "A", "B"]          # sorted by score, zero score dropped
    assert list(net.links.TG) == ["x", "y", "z", "x"]
    # links of TF A form one module (3 nodes, spanning tree with 2 edges); B->x is a singleton
    assert net.modules.module.nunique() == 1 and len(net.modules) == 3
    assert len(net.edges) == 2
    # maximum spanning tree keeps the two strongest of the three co-regulation edges
    w = {(min(u, v), max(u, v)) for u, v in net.edges[["u", "v"]].to_numpy()}
    assert w == {(0, 1), (1, 2)}
    np.testing.assert_allclose(net.adj[0, 1], 0.8 / (1 + 0.4))


def test_coregulated_pairs_jaccard():
    net = {("A", "x"), ("A", "y"), ("B", "x"), ("B", "y"), ("C", "z"), ("A", "z")}
    pairs = coregulated_pairs(net)
    assert ("x", "y") in pairs and ("y", "x") in pairs        # identical regulators
    assert ("x", "z") not in pairs                             # {A,B} vs {A,C}: 1/3
