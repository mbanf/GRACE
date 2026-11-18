"""GRACE MCP server: explore, rebuild, retrain and validate GRACE networks from an MCP client.

Usage:
    python -m grace_mcp --data-dir /path/to/datasets        # or set GRACE_DATA_DIR
The data folder holds the unpacked athaliana_datasets/ and drosophila_datasets/.

Networks are kept by name per species: "published" is the network of the paper; grace_build_network
and grace_train add more (e.g. "published@0.7", "retrained"), which all tools can then query.
"""

from __future__ import annotations

import argparse
import os
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
try:                                           # MCP Python SDK 2.x
    from mcp.server.mcpserver import MCPServer as FastMCP
except ImportError:                            # MCP Python SDK 1.x
    from mcp.server.fastmcp import FastMCP

from grace.pipeline import Prepared, ensemble_network, prepare
from grace.reference import SPECIES, SPECIES_NAMES, published_models, published_network
from grace.train import Model, default_workers, ensemble, train_models
from grace.validation import evidence_sets, validate

mcp = FastMCP("grace", instructions=(
    "GRACE (Banf & Rhee, Sci Rep 2017) refines an initial gene regulatory network (random forest scores "
    "+ conserved motifs) with Markov random fields over co-regulation modules, trained semi-supervised "
    "on a regulatory gold standard and co-functional gene pairs; the final network is the ensemble of "
    "100 bootstrap models. Species: 'athaliana' (AGI ids or gene symbols, e.g. FLC) and 'drosophila' "
    "(FlyBase FBgn ids). Start with grace_summary."))


class State:
    def __init__(self, data_dir: Path | None):
        self.data_dir = data_dir
        self.prepared: dict[str, Prepared] = {}
        self.networks: dict[str, dict[str, pd.DataFrame]] = {s: {} for s in SPECIES}
        self.models: dict[str, dict[str, list[Model]]] = {s: {"published": published_models(s)} for s in SPECIES}

    def get(self, species: str) -> Prepared:
        if species not in SPECIES:
            raise ValueError(f"unknown species {species!r}; use one of {SPECIES}")
        if species not in self.prepared:
            if not self.data_dir or not (self.data_dir / f"{species}_datasets").exists():
                raise RuntimeError(f"{species}_datasets/ not found; start the server with --data-dir or GRACE_DATA_DIR")
            p = prepare(species, self.data_dir)
            self.prepared[species] = p
            self.networks[species]["published"] = ensemble_network(p, self.models[species]["published"])
        return self.prepared[species]

    def network(self, species: str, name: str) -> pd.DataFrame:
        self.get(species)
        if name not in self.networks[species]:
            raise ValueError(f"no network {name!r} for {species}; available: {sorted(self.networks[species])}")
        return self.networks[species][name]


STATE = State(None)


# ---------------------------------------------------------------- gene name handling
def _species_of(gene: str) -> str | None:
    g = gene.upper()
    if g.startswith("FBGN"):
        return "drosophila"
    if len(g) == 9 and g.startswith("AT") and g[2] in "12345CM" and g[3] == "G":
        return "athaliana"
    return None


def _resolve(gene: str, species: str | None) -> tuple[str, str]:
    """(species, gene id) for an id or an A. thaliana gene symbol."""
    sp = species or _species_of(gene)
    if sp == "drosophila":
        return sp, gene if gene.startswith("FBgn") else "FBgn" + gene[4:]
    if sp == "athaliana" and _species_of(gene) == "athaliana":
        return sp, gene.upper()
    # symbol lookup (A. thaliana)
    ds = STATE.get("athaliana").dataset
    hits = [k for k, v in ds.symbols.items() if v.upper() == gene.upper()]
    if not hits:
        raise ValueError(f"unknown gene {gene!r}: use an AGI id (e.g. AT5G10140), a gene symbol (e.g. FLC) or an FBgn id")
    return "athaliana", hits[0]


def _label(species: str, gene: str) -> dict:
    ds = STATE.get(species).dataset
    d = {"id": gene}
    if gene in ds.symbols:
        d["symbol"] = ds.symbols[gene]
    if gene in ds.families:
        d["tf_family"] = ds.families[gene]
    return d


def _overlap(a: set, b: set) -> dict:
    return {"size_a": len(a), "size_b": len(b), "shared": len(a & b), "only_a": len(a - b), "only_b": len(b - a),
            "jaccard": round(len(a & b) / len(a | b), 4) if a | b else 1.0}


# ---------------------------------------------------------------- tools
@mcp.tool()
def grace_summary(species: str | None = None) -> dict:
    """Overview of the GRACE inputs and networks for one species ('athaliana', 'drosophila') or both."""
    out = {}
    for s in [species] if species else SPECIES:
        p = STATE.get(s)
        t = p.training
        out[s] = {
            "species": SPECIES_NAMES[s],
            "candidate_links": len(p.network.links),
            "transcription_factors": len(p.network.tfs),
            "target_genes": len(p.network.tgs),
            "coregulation_modules": int(p.network.modules.module.nunique()),
            "gold_standard_links": len(p.gs_pairs),
            "cofunctional_evidence_pairs": int(t.bp.nnz),
            "networks": {k: len(v) for k, v in STATE.networks[s].items()},
            "model_sets": {k: len(v) for k, v in STATE.models[s].items()},
            "validation_sets": list(p.dataset.validation),
        }
    return out


@mcp.tool()
def grace_query_tf(tf: str, species: str | None = None, network: str = "published", limit: int = 50) -> dict:
    """Targets of a transcription factor in a GRACE network, with ensemble support and initial score."""
    species, tf = _resolve(tf, species)
    p = STATE.get(species)
    net = STATE.network(species, network)
    L = p.network.links
    cand = L[L.TF == tf]
    tg = net[net.TF == tf].sort_values(["support", "score"], ascending=False)
    gs = {g for t, g in p.gs_pairs if t == tf}
    return {"tf": _label(species, tf), "species": species, "network": network,
            "candidate_links": len(cand), "targets_in_network": len(tg),
            "gold_standard_targets": sorted(gs), "gold_standard_recovered": sorted(gs & set(tg.TG)),
            "targets": [{**_label(species, r.TG), "support": round(float(r.support), 3), "score": float(r.score)}
                        for r in tg.head(limit).itertuples()]}


@mcp.tool()
def grace_query_gene(gene: str, species: str | None = None, network: str = "published") -> dict:
    """Regulators of a target gene in a GRACE network."""
    species, gene = _resolve(gene, species)
    p = STATE.get(species)
    net = STATE.network(species, network)
    regs = net[net.TG == gene].sort_values(["support", "score"], ascending=False)
    cand = p.network.links[p.network.links.TG == gene]
    return {"gene": _label(species, gene), "species": species, "network": network,
            "candidate_regulators": len(cand), "regulators_in_network": len(regs),
            "regulators": [{**_label(species, r.TF), "support": round(float(r.support), 3), "score": float(r.score)}
                           for r in regs.itertuples()]}


@mcp.tool()
def grace_query_link(tf: str, target: str, species: str | None = None) -> dict:
    """Everything GRACE knows about one TF -> target link: initial score, co-regulation module, Markov
    random field neighbours, support in each network, and the evidence sets that contain it."""
    species, tf = _resolve(tf, species)
    _, target = _resolve(target, species)
    p = STATE.get(species)
    L = p.network.links
    idx = np.nonzero((L.TF.to_numpy() == tf) & (L.TG.to_numpy() == target))[0]
    if not len(idx):
        return {"tf": _label(species, tf), "target": _label(species, target), "candidate": False,
                "note": "not a candidate link: no initial GRN score above the pre-selection threshold"}
    i = int(idx[0])
    mod = p.network.modules
    row = mod[mod.link == i]
    neighbours = []
    if len(row):
        e = p.network.edges
        nb = e[(e.u == i) | (e.v == i)]
        for r in nb.itertuples():
            j = r.v if r.u == i else r.u
            neighbours.append({"tf": L.TF[j], "target": _label(species, L.TG[j])["id"], "weight": round(float(r.weight), 4)})
    support = {name: (float(df.loc[(df.TF == tf) & (df.TG == target), "support"].iloc[0])
                      if ((df.TF == tf) & (df.TG == target)).any() else 0.0)
               for name, df in STATE.networks[species].items()}
    ev = []
    if (tf, target) in p.dataset.regulatory_evidence:
        ev.append("regulatory gold standard")
    for k, v in p.dataset.validation.items():
        if isinstance(v, set) and (tf, target) in v:
            ev.append(k)
    return {"tf": _label(species, tf), "target": _label(species, target), "candidate": True,
            "initial_score": float(L.v_grn[i]),
            "module": None if not len(row) else {"id": int(row.module.iloc[0]),
                                                  "size": int((mod.module == row.module.iloc[0]).sum())},
            "singleton": not len(row), "spanning_tree_neighbours": neighbours,
            "support_by_network": support, "in_evidence": ev}


@mcp.tool()
def grace_build_network(species: str, threshold: float = 0.5, models: str = "published", name: str | None = None) -> dict:
    """Build an ensemble network from a model set ('published' or a retrained set) keeping links predicted
    by more than `threshold` of the models; stored under `name` (default '<models>@<threshold>')."""
    p = STATE.get(species)
    if models not in STATE.models[species]:
        raise ValueError(f"no model set {models!r}; available: {sorted(STATE.models[species])}")
    net = ensemble_network(p, STATE.models[species][models], threshold=threshold)
    name = name or f"{models}@{threshold:g}"
    STATE.networks[species][name] = net
    got = set(zip(net.TF, net.TG))
    return {"network": name, "links": len(net), "tfs": int(net.TF.nunique()), "targets": int(net.TG.nunique()),
            "vs_published": _overlap(got, published_network(species))}


@mcp.tool()
def grace_train(species: str, n_models: int = 100, seed_offset: int = 0, workers: int | None = None,
                name: str = "retrained") -> dict:
    """Retrain GRACE from scratch (bootstrap + grid search + simulated annealing per model) and build the
    ensemble network. Takes ~20 s (athaliana) or ~90 s (drosophila) for 100 models on 12 cores.
    seed_offset=0 uses the seeds of the original R tutorials."""
    p = STATE.get(species)
    workers = workers or default_workers()
    seeds = [1234 + 158 * i + seed_offset for i in range(1, n_models + 1)]
    t = time.perf_counter()
    models = train_models(p.training, n_workers=workers, seeds=seeds)
    dt = time.perf_counter() - t
    STATE.models[species][name] = models
    net = ensemble_network(p, models)
    STATE.networks[species][name] = net
    df = pd.DataFrame([asdict(m) for m in models])
    ok = df[df.gamma > 0]
    pub = pd.DataFrame([asdict(m) for m in STATE.models[species]["published"]])
    return {"model_set": name, "network": name, "models": len(models), "seconds": round(dt, 1), "workers": workers,
            "median_gamma": float(ok.gamma.median()), "median_lambda": float(ok.lam.median()),
            "published_median_gamma": float(pub.gamma.median()), "published_median_lambda": float(pub.lam.median()),
            "links": len(net), "vs_published": _overlap(set(zip(net.TF, net.TG)), published_network(species))}


@mcp.tool()
def grace_validate(species: str, network: str = "published") -> dict:
    """Enrichment of a network for regulatory and co-functional evidence: TF->target links against
    gold standards / ChIP / AtRegNet, co-regulated target pairs (Jaccard >= 0.5) against co-functional
    pairs, co-localization, Hi-C contacts. Fold change over background and hypergeometric p."""
    p = STATE.get(species)
    net = STATE.network(species, network)
    links, pairs = evidence_sets(p.dataset)
    res = validate(set(zip(net.TF, net.TG)), p.network.tfs, p.network.tgs, links, pairs)
    return {"network": network, "links": len(net),
            "enrichment": [{k: (round(v, 4) if isinstance(v, float) and k == "fold_change" else v)
                            for k, v in asdict(r).items()} for r in res]}


@mcp.tool()
def grace_compare(species: str, network_a: str = "published", network_b: str = "retrained") -> dict:
    """Compare two networks: shared and unique links, and the TFs that gained or lost most targets."""
    a, b = STATE.network(species, network_a), STATE.network(species, network_b)
    sa, sb = set(zip(a.TF, a.TG)), set(zip(b.TF, b.TG))
    ca, cb = a.TF.value_counts(), b.TF.value_counts()
    diff = (cb.reindex(ca.index.union(cb.index), fill_value=0) - ca.reindex(ca.index.union(cb.index), fill_value=0))
    diff = diff[diff != 0].sort_values(key=abs, ascending=False).head(10)
    return {"a": network_a, "b": network_b, **_overlap(sa, sb),
            "tf_target_count_changes": [{**_label(species, t), "change": int(d)} for t, d in diff.items()]}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="grace-mcp", description="GRACE MCP server (stdio)")
    ap.add_argument("--data-dir", default=os.environ.get("GRACE_DATA_DIR"),
                    help="folder with athaliana_datasets/ and drosophila_datasets/ (or GRACE_DATA_DIR)")
    args = ap.parse_args(argv)
    STATE.data_dir = Path(args.data_dir).expanduser() if args.data_dir else None
    mcp.run()


if __name__ == "__main__":
    main()
