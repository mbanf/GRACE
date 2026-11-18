"""Meta co-regulation network and module construction (port of contruct_grace_meta_connectivity_network,
build_link_ranking_connectivity_matrix and prepare_meta_modules in R/GRACE.R).

Nodes are candidate TF->TG links. Two links of the same TF are connected when their targets are
co-functional; the weight is cofunctional(tg1, tg2) / (1 + |v_grn1 - v_grn2|). Each connected component
becomes one module, reduced to its maximum spanning tree.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree

from .mrf import Forest, build_forest


@dataclass
class GraceNetwork:
    links: pd.DataFrame      # TF, TG, v_grn, tf_id, tg_id (0-based); row i = link index i
    tfs: list[str]
    tgs: list[str]
    adj: sp.csr_matrix       # link x link co-regulation weights (symmetric)
    modules: pd.DataFrame    # module, pos (1-based), link
    edges: pd.DataFrame      # module, u, v (link indices), weight
    forest: Forest


def build_links(grn: pd.DataFrame, th_preselection: float = 0.0) -> pd.DataFrame:
    """grn: long table (TF, TG, value) in column-major order of the TF x TG matrix (as R's as.table)."""
    df = grn[grn["value"] > 0].rename(columns={"value": "v_grn"}).reset_index(drop=True)
    # R: order(-v.grn) is a stable radix sort
    df = df.iloc[np.argsort(-df["v_grn"].to_numpy(), kind="stable")].reset_index(drop=True)
    th = max(th_preselection, df["v_grn"].iloc[-1])
    df = df[df["v_grn"] >= th].reset_index(drop=True)
    df["tf_id"] = pd.factorize(df["TF"])[0]
    df["tg_id"] = pd.factorize(df["TG"])[0]
    return df[["TF", "TG", "v_grn", "tf_id", "tg_id"]]


def build_coregulation(links: pd.DataFrame, cofun: sp.csr_matrix) -> sp.csr_matrix:
    """cofun: TG x TG co-functional matrix indexed by tg_id."""
    n = len(links)
    v = links["v_grn"].to_numpy()
    tg = links["tg_id"].to_numpy()
    rows, cols, vals = [], [], []
    for _, idx in links.groupby("tf_id", sort=False).indices.items():
        if len(idx) < 2:
            continue
        idx = np.sort(idx)
        a, b = np.triu_indices(len(idx), k=1)
        i, j = idx[a], idx[b]
        rn = np.asarray(cofun[tg[i], tg[j]]).ravel()
        keep = rn != 0
        if keep.any():
            rows.append(i[keep]); cols.append(j[keep])
            vals.append(rn[keep] / (1.0 + np.abs(v[i[keep]] - v[j[keep]])))
    if not rows:
        return sp.csr_matrix((n, n))
    r, c, w = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
    m = sp.coo_matrix((w, (r, c)), shape=(n, n))
    return (m + m.T).tocsr()


def build_modules(adj: sp.csr_matrix):
    """Connected components of non-singleton links, each reduced to its maximum spanning tree."""
    n = adj.shape[0]
    deg = np.diff(adj.indptr)
    nonsingle = np.nonzero(deg > 0)[0]
    sub = adj[nonsingle][:, nonsingle]
    _, lab = connected_components(sub, directed=False)
    # number components by their lowest link index (igraph order)
    first = {}
    for k, l in enumerate(lab):
        first.setdefault(l, k)
    remap = {l: r for r, l in enumerate(sorted(first, key=first.get))}
    mod = np.array([remap[l] for l in lab])

    # maximum spanning tree = minimum spanning tree on (1 - w); all w in (0, 1) so 1 - w > 0
    inv = sub.copy()
    inv.data = 1.0 - inv.data
    mst = minimum_spanning_tree(inv).tocoo()
    u, w_ = nonsingle[mst.row], nonsingle[mst.col]
    weight = np.asarray(adj[u, w_]).ravel()
    edges = pd.DataFrame({"module": mod[mst.row], "u": u, "v": w_, "weight": weight})

    modules = pd.DataFrame({"module": mod, "link": nonsingle})
    modules = modules.sort_values(["module", "link"], kind="stable").reset_index(drop=True)
    modules["pos"] = modules.groupby("module").cumcount() + 1
    return modules, edges


def forest_from(modules: pd.DataFrame, edges: pd.DataFrame, n_links: int) -> Forest:
    pos = dict(zip(modules["link"], modules["pos"]))
    return build_forest(n_links, modules["module"].to_numpy(), modules["link"].to_numpy(), modules["pos"].to_numpy(),
                        edges["module"].to_numpy(), edges["u"].map(pos).to_numpy(), edges["v"].map(pos).to_numpy(),
                        edges["weight"].to_numpy())


def build_network(grn: pd.DataFrame, cofun: pd.DataFrame, th_preselection: float = 0.0) -> GraceNetwork:
    """grn: long (TF, TG, value) table; cofun: long (row, col, value) TG x TG co-functional table."""
    links = build_links(grn, th_preselection)
    tfs = links.drop_duplicates("tf_id")["TF"].tolist()
    tgs = links.drop_duplicates("tg_id")["TG"].tolist()
    tg_index = {g: k for k, g in enumerate(tgs)}
    cf = cofun[cofun["row"].isin(tg_index) & cofun["col"].isin(tg_index)]
    m = sp.coo_matrix((cf["value"].to_numpy(), (cf["row"].map(tg_index), cf["col"].map(tg_index))),
                      shape=(len(tgs), len(tgs))).tocsr()
    adj = build_coregulation(links, m)
    modules, edges = build_modules(adj)
    return GraceNetwork(links, tfs, tgs, adj, modules, edges, forest_from(modules, edges, len(links)))
