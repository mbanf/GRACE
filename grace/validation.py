"""Enrichment of a GRACE network for independent evidence (as in the tutorials' validation).

Regulatory link sets (gold standard, ChIP, AtRegNet) are compared with the predicted TF->TG links;
gene-pair sets (co-functional pairs, co-localization, Hi-C contacts, metabolic pathways) are compared with
co-regulated target pairs, i.e. targets whose regulator sets have Jaccard similarity >= 0.5.
Fold change = observed rate / background rate; p = one-sided hypergeometric test. Backgrounds are all
candidate TF x TG pairs (links) or all target x target pairs (gene pairs), as in the R tutorials.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.stats import hypergeom

from .datasets import GroupPairs

Pairs = set[tuple[str, str]]


@dataclass
class Enrichment:
    evidence: str
    kind: str          # "links" or "gene pairs"
    n_predicted: int   # predicted links or co-regulated pairs
    hits: int
    in_background: int
    background: int
    fold_change: float
    p_value: float


def coregulated_pairs(network: Pairs) -> Pairs:
    """Ordered target pairs (a != b) whose regulator sets have Jaccard >= 0.5."""
    if not network:
        return set()
    tfs = sorted({t for t, _ in network})
    tgs = sorted({g for _, g in network})
    ti = {t: k for k, t in enumerate(tfs)}
    gi = {g: k for k, g in enumerate(tgs)}
    r, c = zip(*[(gi[g], ti[t]) for t, g in network])
    M = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(tgs), len(tfs)))
    A = (M @ M.T).tocoo()
    b = np.asarray(M.sum(axis=1)).ravel()
    off = A.row != A.col
    i, j, a = A.row[off], A.col[off], A.data[off]
    keep = a / (b[i] + b[j] - a) >= 0.5
    return {(tgs[x], tgs[y]) for x, y in zip(i[keep], j[keep])}


def _result(name, kind, n, hits, pop, background) -> Enrichment:
    fold = (hits / n) / (pop / background) if n and pop else float("nan")
    p = float(hypergeom.sf(hits - 1, background, pop, n)) if n and pop else float("nan")
    return Enrichment(name, kind, n, hits, pop, background, fold, p)


def _pairs_test(name, kind, predicted: Pairs, evidence, universe_a, universe_b) -> Enrichment:
    ua, ub = set(universe_a), set(universe_b)
    background = len(ua) * len(ub)
    if isinstance(evidence, GroupPairs):            # gene pairs sharing a compartment / tissue
        pop = evidence.count_within(ua)
        hits = evidence.hits(predicted)
    else:
        pop_set = {p for p in evidence if p[0] in ua and p[1] in ub and p[0] != p[1]}
        pop, hits = len(pop_set), len(predicted & pop_set)
    return _result(name, kind, len(predicted), hits, pop, background)


def validate(network: Pairs, tfs, tgs, link_sets: dict, pair_sets: dict) -> list[Enrichment]:
    co = coregulated_pairs(network)
    out = [_pairs_test(k, "links", network, v, tfs, tgs) for k, v in link_sets.items()]
    out += [_pairs_test(k, "gene pairs", co, v, tgs, tgs) for k, v in pair_sets.items()]
    return out


def evidence_sets(dataset) -> tuple[dict, dict]:
    """Split a Dataset's evidence into link sets and gene-pair sets."""
    links = {"regulatory gold standard (training)": dataset.regulatory_evidence}
    pairs = {"co-functional pairs (training)": dataset.cofunctional_evidence}
    for k, v in dataset.validation.items():
        (links if k in ("AtRegNet", "ChIP") else pairs)[k] = v
    return links, pairs
