"""Markov Random Field re-evaluation of candidate regulatory links (port of GRACE v1 R code).

Each co-regulation module is a maximum spanning tree over candidate links (nodes = TF->TG links).
Node potentials: [exp(v_grn - gamma), 1] for states (link, no link); links with v_grn > alpha are clamped
to "link". Edge potentials are Potts: exp(lambda * w) on agreement, 1 otherwise. Inference is exact
sum-product belief propagation on the trees, done in log-odds space. A link is predicted when its
marginal P(link) > 0.5, i.e. log-odds > 0. This matches CRF::infer.conditional(crf, clamped, infer.tree)
as used in reevaluate_links_with_MarkovRandomFields() (R/GRACE.R).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit


@dataclass
class Forest:
    """All module trees stacked into one forest, nodes indexed by link index (0-based)."""

    order: np.ndarray        # node indices in BFS order (roots first within each module)
    parent: np.ndarray       # parent node per node, -1 for roots (indexed by node)
    parent_w: np.ndarray     # weight of the edge to the parent (indexed by node)
    module_of: np.ndarray    # module index per node, -1 if the node is in no module
    module_root: np.ndarray  # root node per module
    n_modules: int


def build_forest(n_links: int, node_module: np.ndarray, node_link: np.ndarray, node_pos: np.ndarray,
                 edge_module: np.ndarray, edge_from: np.ndarray, edge_to: np.ndarray, edge_w: np.ndarray) -> Forest:
    """Build a rooted forest from module node/edge tables (positions are 1-based within a module)."""
    n_modules = int(node_module.max()) + 1 if len(node_module) else 0
    # map (module, pos) -> link index
    key = {}
    for m, p, l in zip(node_module, node_pos, node_link):
        key[(int(m), int(p))] = int(l)
    adj: dict[int, list[tuple[int, float]]] = {}
    for m, a, b, w in zip(edge_module, edge_from, edge_to, edge_w):
        u, v = key[(int(m), int(a))], key[(int(m), int(b))]
        adj.setdefault(u, []).append((v, float(w)))
        adj.setdefault(v, []).append((u, float(w)))

    module_of = np.full(n_links, -1, dtype=np.int64)
    module_of[node_link] = node_module
    parent = np.full(n_links, -1, dtype=np.int64)
    parent_w = np.zeros(n_links, dtype=np.float64)
    module_root = np.full(n_modules, -1, dtype=np.int64)
    order: list[int] = []
    seen = np.zeros(n_links, dtype=bool)
    # root = first node (pos 1) of each module, matching node order of the R export
    roots = {}
    for m, p, l in zip(node_module, node_pos, node_link):
        if p == 1:
            roots[int(m)] = int(l)
    for m in range(n_modules):
        r = roots[m]
        module_root[m] = r
        seen[r] = True
        queue = [r]
        head = 0
        while head < len(queue):
            u = queue[head]
            head += 1
            order.append(u)
            for v, w in adj.get(u, ()):
                if not seen[v]:
                    seen[v] = True
                    parent[v] = u
                    parent_w[v] = w
                    queue.append(v)
    return Forest(np.asarray(order, dtype=np.int64), parent, parent_w, module_of, module_root, n_modules)


@njit(cache=True)
def _logaddexp(a, b):
    if a > b:
        return a + np.log1p(np.exp(b - a))
    return b + np.log1p(np.exp(a - b))


@njit(cache=True)
def _msg(h, J):
    """Log-odds message through a Potts edge with coupling J given sender log-odds h."""
    return _logaddexp(J + h, 0.0) - _logaddexp(h, J)


@njit(cache=True)
def _infer_forest(order, parent, parent_w, active_node, v_grn, alpha, gamma, lam, out):
    """Exact BP over all active nodes in `order`. Writes 0/1 predictions into out[node]."""
    n = v_grn.shape[0]
    up = np.zeros(n)        # unary + messages from children
    m_up = np.zeros(n)      # message node -> parent
    total = np.zeros(n)     # final log-odds
    clamped = np.zeros(n, dtype=np.bool_)
    k = order.shape[0]
    for t in range(k):
        u = order[t]
        if active_node[u]:
            up[u] = v_grn[u] - gamma
            clamped[u] = v_grn[u] > alpha
    # upward pass
    for t in range(k - 1, -1, -1):
        u = order[t]
        if not active_node[u]:
            continue
        p = parent[u]
        if p >= 0:
            J = lam * parent_w[u]
            if clamped[u]:
                m = J
            else:
                m = _msg(up[u], J)
            m_up[u] = m
            up[p] += m
    # downward pass
    for t in range(k):
        u = order[t]
        if not active_node[u]:
            continue
        p = parent[u]
        if p < 0:
            total[u] = up[u]
        else:
            J = lam * parent_w[u]
            if clamped[p]:
                md = J
            else:
                md = _msg(total[p] - m_up[u], J)
            total[u] = up[u] + md
        if clamped[u] or total[u] > 0.0:
            out[u] = 1.0
        else:
            out[u] = 0.0


def predict_links(forest: Forest, v_grn: np.ndarray, adj_indptr: np.ndarray, adj_indices: np.ndarray,
                  alpha: float, gamma: float, lam: float,
                  active_links: np.ndarray | None = None, active_modules: np.ndarray | None = None) -> np.ndarray:
    """Per-link prediction values (0, 1, or 2), as in reevaluate_links_with_MarkovRandomFields(b.prim=TRUE).

    active_links: bool mask of links kept in a bootstrap sample (default all). Singletons are links with no
                  co-regulation partner among the active links; they are predicted when exp(v - gamma) > 1.
    active_modules: bool mask of modules kept (default all). Module predictions use the full module tree.
    """
    n = v_grn.shape[0]
    if active_links is None:
        active_links = np.ones(n, dtype=bool)
    if active_modules is None:
        active_modules = np.ones(forest.n_modules, dtype=bool)
    out = np.zeros(n)
    node_active = np.zeros(n, dtype=np.bool_)
    in_mod = forest.module_of >= 0
    node_active[in_mod] = active_modules[forest.module_of[in_mod]]
    _infer_forest(forest.order, forest.parent, forest.parent_w, node_active, v_grn, alpha, gamma, lam, out)
    # singletons: active links whose active co-regulation degree is 0
    single = _singletons(adj_indptr, adj_indices, active_links)
    s = single & (np.exp(v_grn - gamma) > 1.0)
    out[s] += 1.0
    return out


@njit(cache=True)
def _singletons_nb(indptr, indices, active):
    n = active.shape[0]
    res = np.zeros(n, dtype=np.bool_)
    for i in range(n):
        if not active[i]:
            continue
        deg = 0
        for k in range(indptr[i], indptr[i + 1]):
            if active[indices[k]]:
                deg += 1
                break
        res[i] = deg == 0
    return res


def _singletons(indptr, indices, active):
    return _singletons_nb(indptr, indices, np.ascontiguousarray(active, dtype=np.bool_))
