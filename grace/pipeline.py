"""End-to-end GRACE v1 pipeline from the tutorial datasets (no R needed).

Mirrors GRACE_tutorial_athaliana.R / GRACE_tutorial_drosophila.R up to model training: hyperparameter grids,
pre-selection, meta network and modules, regulatory gold standard (truereg) and co-functional evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from .datasets import Dataset, load_dataset
from .network import GraceNetwork, build_network
from .train import TrainingData

# settings of the two R tutorials
SETTINGS = {
    "athaliana": dict(max_call=250, coreg_column="n.coreg.tg_pairs.BP"),
    "drosophila": dict(max_call=200, coreg_column="n.coreg.bp_tg_pairs"),
}
N_SETS = 20
LAMBDA_GRID = np.array([0.01, 0.5, 1.0, 1.5, 2.0, 2.5])
N_SAMPLE_SIZE = 0.632
TH_PERCENTAGE = 0.99
BETA = 1.0


@dataclass
class Prepared:
    dataset: Dataset
    network: GraceNetwork
    training: TrainingData
    gs_pairs: set           # truereg pairs (TF, TG)
    constants: dict


def gamma_grid(grn_values: np.ndarray, n_sets: int = N_SETS) -> np.ndarray:
    pos = grn_values[grn_values > 0]
    upper = float(np.quantile(pos, 0.99999999, method="linear"))   # R quantile type 7
    lower = float(pos.min())
    return np.array([upper - (upper - lower) / n_sets * m for m in range(1, n_sets + 1)])


def normalization_constants(rate: "pd.DataFrame", coreg_column: str) -> tuple[float, float]:
    """ratio.gs.grn_vs_bp.normalization and max.coreg.bp from the precomputed rate-density curve."""
    precision = rate["n.gs.grn"].to_numpy() / rate["n.pred"].to_numpy()
    return 1.0 / float(np.max(precision)), float(rate[coreg_column].max())


def prepare(species: str, data_root: str | Path) -> Prepared:
    ds = load_dataset(species, Path(data_root))
    cfg = SETTINGS[species]
    v = ds.grn["value"].to_numpy(float)
    th_preselection = float(v[v > 0].min())
    net = build_network(ds.grn, ds.cofunctional_network, th_preselection=th_preselection)

    L = net.links
    gs_pairs = {(a, b) for a, b in zip(L.TF, L.TG) if (a, b) in ds.regulatory_evidence}
    gs_link = np.array([1.0 if (a, b) in ds.regulatory_evidence else 0.0 for a, b in zip(L.TF, L.TG)])
    gs_tfs = np.zeros(len(net.tfs), dtype=bool)
    gs_tfs[L.tf_id.to_numpy()[gs_link > 0]] = True

    tg_index = {g: k for k, g in enumerate(net.tgs)}
    bp_pairs = [(tg_index[a], tg_index[b]) for a, b in ds.cofunctional_evidence if a in tg_index and b in tg_index]
    r, c = (np.array(x, dtype=int) for x in zip(*bp_pairs)) if bp_pairs else (np.array([], int), np.array([], int))
    bp = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(net.tgs), len(net.tgs)))

    ratio, max_coreg = normalization_constants(ds.rate_density, cfg["coreg_column"])
    first = net.modules[net.modules["pos"] == 1].sort_values("module")["link"].to_numpy()
    td = TrainingData(
        tf_id=L.tf_id.to_numpy(), tg_id=L.tg_id.to_numpy(), v_grn=L.v_grn.to_numpy(),
        n_tfs=len(net.tfs), n_tgs=len(net.tgs), adj_indptr=net.adj.indptr, adj_indices=net.adj.indices,
        forest=net.forest, module_first_tf=L.tf_id.to_numpy()[first], gs_link=gs_link, gs_tfs=gs_tfs, bp=bp,
        ratio_norm=ratio, max_coreg_bp=max_coreg, gamma_grid=gamma_grid(v), lambda_grid=LAMBDA_GRID,
        n_sample_size=N_SAMPLE_SIZE, max_call=cfg["max_call"], beta=BETA, th_percentage=TH_PERCENTAGE)
    constants = dict(ratio_gs_grn_vs_bp_normalization=ratio, max_coreg_bp=max_coreg, th_preselection=th_preselection,
                     max_call=cfg["max_call"])
    return Prepared(ds, net, td, gs_pairs, constants)


def ensemble_network(p: Prepared, models, threshold: float = 0.5):
    """Ensemble network as a table (TF, TG, score, support): links predicted by more than `threshold`
    of the models. support = fraction of models predicting the link; score = initial GRN score."""
    import pandas as pd

    from .train import ensemble

    frac = ensemble(p.training, models)
    keep = frac > threshold
    L = p.network.links
    df = pd.DataFrame({"TF": L.TF[keep].to_numpy(), "TG": L.TG[keep].to_numpy(),
                       "score": L.v_grn[keep].to_numpy(), "support": frac[keep]})
    return df.sort_values(["support", "score"], ascending=False, kind="stable").reset_index(drop=True)
