"""Semi-supervised model learning (port of evaluate_mode_parameter, gridSearchOptimization,
generate_bootstrapping_sample and compute_model_assessment_parallel_optimized).

One model = one bootstrap sample of TFs + regulatory gold standard, a grid search over (gamma, lambda),
then generalized simulated annealing (scipy.optimize.dual_annealing, a port of R's GenSA) in the
neighbourhood of the best grid cell. The objective is the f-score combining normalized precision on the
regulatory gold standard with recovered co-functional target pairs.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.optimize import dual_annealing

from .mrf import Forest, predict_links


@dataclass
class TrainingData:
    tf_id: np.ndarray            # per link
    tg_id: np.ndarray            # per link
    v_grn: np.ndarray            # per link
    n_tfs: int
    n_tgs: int
    adj_indptr: np.ndarray
    adj_indices: np.ndarray
    forest: Forest
    module_first_tf: np.ndarray  # TF id of the first node of each module (R keeps a module if this TF is sampled)
    gs_link: np.ndarray          # per link: regulatory gold standard value (truereg[TF, TG])
    gs_tfs: np.ndarray           # bool per TF: has at least one gold-standard link
    bp: sp.csr_matrix            # TG x TG co-functional evidence (indexed by tg_id; 0 for TGs not in evidence)
    ratio_norm: float            # ratio.gs.grn_vs_bp.normalization
    max_coreg_bp: float
    gamma_grid: np.ndarray       # descending
    lambda_grid: np.ndarray
    n_sample_size: float = 0.632
    max_call: int = 200
    beta: float = 1.0
    th_percentage: float = 0.99
    normalize_precision: bool = True


@dataclass
class Model:
    alpha: float
    gamma: float
    lam: float
    fscore: float
    n_pred: float
    n_gs: float
    n_coreg_bp: float
    seed: int
    grid_best: tuple = field(default=())


# ---------------------------------------------------------------- assessment
def assess(d: TrainingData, pred: np.ndarray, gs_link: np.ndarray, bp: sp.csr_matrix):
    """n.pred, n.gs.grn, n.coreg.bp_tg_pairs as in compute_model_assessment_parallel_optimized."""
    nz = np.nonzero(pred)[0]
    n_pred = float(pred[nz].sum())
    if n_pred == 0:
        return 0.0, 0.0, 0.0
    n_gs = float((pred[nz] * gs_link[nz]).sum())
    # targets x TFs prediction matrix, Jaccard similarity between targets over their regulators
    M = sp.csr_matrix((pred[nz], (d.tg_id[nz], d.tf_id[nz])), shape=(d.n_tgs, d.n_tfs))
    A = (M @ M.T).tocoo()
    b = np.asarray(M.sum(axis=1)).ravel()
    off = A.row != A.col
    i, j, a = A.row[off], A.col[off], A.data[off]
    keep = a > 0
    i, j, a = i[keep], j[keep], a[keep]
    J = a / (b[i] + b[j] - a)
    sel = J >= 0.5
    n_coreg_bp = float(np.asarray(bp[i[sel], j[sel]]).sum()) if sel.any() else 0.0
    return n_pred, n_gs, n_coreg_bp


def fscore(d: TrainingData, n_pred, n_gs, n_coreg_bp):
    y1 = n_gs / n_pred if n_pred > 0 else 0.0
    if d.normalize_precision:
        y1 *= d.ratio_norm
    y2 = n_coreg_bp / d.max_coreg_bp if d.max_coreg_bp > 0 else 0.0
    den = d.beta ** 2 * y1 + y2
    return (1 + d.beta ** 2) * y1 * y2 / den if den > 0 else 0.0


# ---------------------------------------------------------------- bootstrap
@dataclass
class Bootstrap:
    active_links: np.ndarray
    active_modules: np.ndarray
    gs_link: np.ndarray
    bp: sp.csr_matrix


def bootstrap_sample(d: TrainingData, rng: np.random.Generator) -> Bootstrap:
    gs_idx = np.nonzero(d.gs_link > 0)[0]                      # gold-standard links (rows of df.GS)
    n_gs = len(gs_idx)
    sampled = rng.choice(n_gs, size=int(round(n_gs * d.n_sample_size)), replace=False)
    tfs_gs_sample = np.unique(d.tf_id[gs_idx[sampled]])
    n_extra = int(round(d.n_tfs * d.n_sample_size)) - len(tfs_gs_sample)
    non_gs_tfs = np.nonzero(~d.gs_tfs)[0]
    extra = rng.choice(non_gs_tfs, size=max(0, min(n_extra, len(non_gs_tfs))), replace=False)
    tf_keep = np.zeros(d.n_tfs, dtype=bool)
    tf_keep[tfs_gs_sample] = True
    tf_keep[extra] = True

    active_links = tf_keep[d.tf_id]
    active_modules = tf_keep[d.module_first_tf]
    gs_link = np.where(active_links, d.gs_link, 0.0)
    # co-functional evidence: zero the block of targets not covered by the sample
    tg_in = np.zeros(d.n_tgs, dtype=bool)
    tg_in[d.tg_id[active_links]] = True
    bp = d.bp.tocoo()
    keep = tg_in[bp.row] | tg_in[bp.col]
    bp = sp.csr_matrix((bp.data[keep], (bp.row[keep], bp.col[keep])), shape=d.bp.shape)
    return Bootstrap(active_links, active_modules, gs_link, bp)


# ---------------------------------------------------------------- objective
def _objective(d: TrainingData, bs: Bootstrap, gamma: float, lam: float):
    pred = predict_links(d.forest, d.v_grn, d.adj_indptr, d.adj_indices, gamma, gamma, lam,
                         active_links=bs.active_links, active_modules=bs.active_modules)
    if pred.sum() <= 0:
        return -100.0, (0.0, 0.0, 0.0)
    counts = assess(d, pred, bs.gs_link, bs.bp)
    return fscore(d, *counts), counts


def _r_quantile7(x, q):
    """R quantile(type = 7)."""
    return float(np.quantile(x, q, method="linear"))


def grid_search(d: TrainingData, bs: Bootstrap):
    G, L = np.meshgrid(d.gamma_grid, d.lambda_grid, indexing="xy")   # expand.grid: gamma varies fastest
    gammas, lambdas = G.ravel(), L.ravel()
    f = np.array([_objective(d, bs, g, l)[0] for g, l in zip(gammas, lambdas)])
    q = _r_quantile7(f, d.th_percentage)
    sel = f >= q
    g_sel, l_sel, f_sel = gammas[sel], lambdas[sel], f[sel]
    # duplicated f-scores: keep only the smallest lambda among them
    keep = np.ones(len(f_sel), dtype=bool)
    for val in np.unique(f_sel):
        same = f_sel == val
        if same.sum() > 1:
            lmin = l_sel[same].min()
            keep &= ~(same & (l_sel != lmin))
    g_sel, l_sel = g_sel[keep], l_sel[keep]
    return _bounds(np.unique(g_sel), d.gamma_grid, descending=True), _bounds(np.unique(l_sel), d.lambda_grid, descending=False)


def _bounds(selection, grid, descending):
    """(start, lower, upper) as in gridSearchOptimization; gamma grid is descending, lambda ascending."""
    if len(selection) == 1:
        v = float(selection[0])
        k = int(np.nonzero(grid == v)[0][0])
        prev_ = grid[k - 1] if k > 0 else v
        next_ = grid[k + 1] if k < len(grid) - 1 else v
        lo, hi = (next_, prev_) if descending else (prev_, next_)
        return v, float(min(lo, v)), float(max(hi, v))
    return float(selection.min()), float(selection.min()), float(selection.max())


def train_model(d: TrainingData, seed: int) -> Model:
    rng = np.random.default_rng(seed)
    bs = bootstrap_sample(d, rng)
    (g0, g_lo, g_hi), (l0, l_lo, l_hi) = grid_search(d, bs)
    eps = 1e-12
    bounds = [(g_lo, g_hi if g_hi > g_lo else g_lo + eps), (l_lo, l_hi if l_hi > l_lo else l_lo + eps)]
    res = dual_annealing(lambda w: -_objective(d, bs, w[0], w[1])[0], bounds, x0=np.array([g0, l0]),
                         maxfun=d.max_call, seed=int(rng.integers(2**31 - 1)))
    gamma, lam = float(res.x[0]), float(res.x[1])
    # final assessment on the full gold standard, with the bootstrap network
    pred = predict_links(d.forest, d.v_grn, d.adj_indptr, d.adj_indices, gamma, gamma, lam,
                         active_links=bs.active_links, active_modules=bs.active_modules)
    if pred.sum() <= 0:
        return Model(-1, -1, -1, -100, -1, -1, -1, seed, (g0, l0))
    counts = assess(d, pred, d.gs_link, d.bp)
    return Model(gamma, gamma, lam, fscore(d, *counts), *counts, seed=seed, grid_best=(g0, l0))


# ---------------------------------------------------------------- many models in parallel
_D: TrainingData | None = None


def _init(d):
    global _D
    _D = d


def _work(seed):
    return train_model(_D, seed)


def default_workers() -> int:
    """Leave two cores free so the machine stays responsive."""
    return max(1, (os.cpu_count() or 2) - 2)


def train_models(d: TrainingData, n_models: int = 100, n_workers: int | None = None,
                 seeds: list[int] | None = None) -> list[Model]:
    seeds = seeds or [1234 + 158 * i for i in range(1, n_models + 1)]   # same seeds as the R tutorials
    n_workers = n_workers or default_workers()
    if n_workers == 1:
        return [train_model(d, s) for s in seeds]
    with ProcessPoolExecutor(n_workers, initializer=_init, initargs=(d,)) as ex:
        return list(ex.map(_work, seeds))


def ensemble(d: TrainingData, models: list[Model], threshold: float = 0.5) -> np.ndarray:
    """Fraction of models predicting each link on the full network; links with fraction > threshold are kept."""
    acc = np.zeros(len(d.v_grn))
    ok = [m for m in models if m.gamma > 0]
    for m in ok:
        acc += predict_links(d.forest, d.v_grn, d.adj_indptr, d.adj_indices, m.alpha, m.gamma, m.lam)
    return acc / len(ok)
