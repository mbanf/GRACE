"""Tutorial dataset loading (port of request_dataset() in tutorials/GRACE_tutorial_datasets.R).

Matrices are kept as long tables or pair sets keyed by gene names, so no dense genome x genome matrices are
built. R quirks that change results are reproduced: read.table treats '#' as a comment, reshape2::acast
counts entries when (row, col) pairs repeat, and AraNet/FlyNet scores are rescaled to [0, 1] by dividing
by the maximum over the TG x TG matrix (its minimum is 0).
"""

from __future__ import annotations

import csv
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

Pairs = set[tuple[str, str]]


# ---------------------------------------------------------------- R-compatible readers
def read_table(path: Path, header: bool = False, comment: str | None = "#") -> pd.DataFrame:
    """R read.table(sep = "\\t", quote = "", stringsAsFactors = FALSE) for tab-separated files."""
    # latin-1 maps bytes 1:1 (some files hold Mac Roman text); gene identifiers are ASCII
    df = pd.read_csv(path, sep="\t", header=0 if header else None, comment=comment, quoting=csv.QUOTE_NONE,
                     dtype=str, keep_default_na=False, na_values=["NA"], encoding="latin-1")
    if not header:
        df.columns = [f"V{i + 1}" for i in range(df.shape[1])]
    return df


def acast_pairs(df: pd.DataFrame, row: str, col: str) -> Pairs:
    """Non-zero pattern of acast(df, row ~ col) followed by mat[mat != 0] <- 1 (presence of a pair)."""
    return set(zip(df[row], df[col]))


def read_rds_matrix(path: Path) -> tuple[np.ndarray, list[str], list[str], bool]:
    """Read a dense matrix or Matrix::dgCMatrix with dimnames from .rds.

    Returns (values, rownames, colnames, is_sparse). For sparse input, values are COO triplets
    (i, j, x) in column-major order; for dense input, the 2D array.
    """
    import rdata

    def dgc(obj, attrs):
        return {"i": np.asarray(attrs["i"]), "p": np.asarray(attrs["p"]), "x": np.asarray(attrs["x"]),
                "Dim": np.asarray(attrs["Dim"]), "Dimnames": attrs["Dimnames"]}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj = rdata.read_rds(path, constructor_dict={**rdata.conversion.DEFAULT_CLASS_MAP, "dgCMatrix": dgc})
    if isinstance(obj, dict) and "p" in obj:
        nrow, ncol = map(int, obj["Dim"])
        rows, cols = list(obj["Dimnames"][0]), list(obj["Dimnames"][1])
        j = np.repeat(np.arange(ncol), np.diff(obj["p"]))
        return np.c_[obj["i"], j, obj["x"]], rows, cols, True
    # dense matrix -> xarray DataArray
    arr = np.asarray(obj)
    rows = [str(v) for v in obj.coords[obj.dims[0]].values]
    cols = [str(v) for v in obj.coords[obj.dims[1]].values]
    return arr, rows, cols, False


def matrix_long(values, rows, cols, sparse) -> pd.DataFrame:
    """Non-zero entries as (TF, TG, value) in column-major order, as R's as.table(matrix)."""
    if sparse:
        i, j, x = values[:, 0].astype(int), values[:, 1].astype(int), values[:, 2]
        keep = x != 0
        i, j, x = i[keep], j[keep], x[keep]
    else:
        jj, ii = np.nonzero(values.T)          # column-major
        i, j, x = ii, jj, values[ii, jj]
    return pd.DataFrame({"TF": np.asarray(rows, dtype=object)[i], "TG": np.asarray(cols, dtype=object)[j], "value": x})


# ---------------------------------------------------------------- dataset container
@dataclass
class Dataset:
    species: str
    tfs: list[str]                       # v.tfs (rows of the GRN matrix)
    tgs: list[str]                       # v.tgs (columns of the GRN matrix)
    grn: pd.DataFrame                    # non-zero GRN scores (TF, TG, value), column-major order
    cofunctional_network: pd.DataFrame   # (row, col, value) TG x TG, rescaled to [0, 1]
    regulatory_evidence: Pairs           # mat.GS (TF, TG) pairs
    cofunctional_evidence: Pairs         # mat.GO.tgs (TG, TG) pairs
    rate_density: pd.DataFrame           # precomputed comparison curve (df.rate_density)
    validation: dict = field(default_factory=dict)          # name -> Pairs or GroupPairs
    symbols: dict[str, str] = field(default_factory=dict)   # gene id -> gene symbol (A. thaliana)
    families: dict[str, str] = field(default_factory=dict)  # TF id -> TF family (A. thaliana)


def _rescale(df: pd.DataFrame, n_cells: int) -> pd.DataFrame:
    """mapping_intervals(mat, min(mat), max(mat), 0, 1) over a matrix where df holds the non-zero cells."""
    mn = min(df["value"].min(), 0.0) if len(df) < n_cells else df["value"].min()
    mx = df["value"].max()
    out = df.copy()
    out["value"] = (df["value"] - mn) / (mx - mn)
    return out


def _restrict(df: pd.DataFrame, a: set, b: set, ca="gn.a", cb="gn.b") -> pd.DataFrame:
    return df[df[ca].isin(a) & df[cb].isin(b)]


class GroupPairs:
    """Gene pairs (a != b) that share at least one group (compartment, tissue), stored as memberships.

    Equivalent to R's jaccard_count(membership) >= 1 with a zero diagonal, without materializing the
    (possibly millions of) pairs.
    """

    def __init__(self, membership: dict[str, set[str]]):
        import scipy.sparse as sp

        self.genes = sorted(g for g, c in membership.items() if c)
        groups = sorted({c for g in self.genes for c in membership[g]})
        self.index = {g: k for k, g in enumerate(self.genes)}
        gi = {c: k for k, c in enumerate(groups)}
        rows, cols = zip(*[(self.index[g], gi[c]) for g in self.genes for c in membership[g]]) if self.genes else ((), ())
        self.B = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(self.genes), len(groups)))

    def _shares(self, a_idx, b_idx) -> np.ndarray:
        return np.asarray(self.B[a_idx].multiply(self.B[b_idx]).sum(axis=1)).ravel() > 0

    def count_within(self, genes) -> int:
        idx = [self.index[g] for g in genes if g in self.index]
        sub = self.B[idx]
        shared = (sub @ sub.T).tocoo()
        return int(((shared.data > 0) & (shared.row != shared.col)).sum())

    def hits(self, pairs) -> int:
        pairs = [(self.index[a], self.index[b]) for a, b in pairs if a != b and a in self.index and b in self.index]
        if not pairs:
            return 0
        a, b = map(np.array, zip(*pairs))
        return int(self._shares(a, b).sum())

    def __len__(self) -> int:
        return self.count_within(self.genes)


def _colocalization(gene_to_compartments: dict[str, set[str]], genes: set[str]) -> GroupPairs:
    """jaccard_count >= 1 with zero diagonal: pairs of genes sharing at least one compartment."""
    return GroupPairs({g: c for g, c in gene_to_compartments.items() if g in genes})


# ---------------------------------------------------------------- A. thaliana
def load_athaliana(root: Path) -> Dataset:
    d = Path(root) / "athaliana_datasets"

    # validation: AtRegNet confirmed links (TFs with < 2000 targets)
    rn = read_table(d / "reg_net.txt")
    rn = rn[rn["V8"] == "Confirmed"]
    rn = pd.DataFrame({"TF": rn["V2"].str.upper(), "Target": rn["V5"].str.upper()})
    atregnet = acast_pairs(rn, "TF", "Target")
    n_targets = pd.Series([t for t, _ in atregnet]).value_counts()
    atregnet = {(t, g) for t, g in atregnet if n_targets[t] < 2000}

    atrm = read_table(d / "ATRM.txt", header=True)
    reg_evidence = acast_pairs(atrm, "TF", "Target")

    values, tfs, tgs, sparse = read_rds_matrix(d / "mat.conserved_motifs_plus_genie3_development.rds")
    grn = matrix_long(values, tfs, tgs, sparse)
    gns = set(tfs) | set(tgs)
    tg_set = set(tgs)

    gs = read_table(d / "AraNet.v2.goldstandard.txt", comment="#")
    gs.columns = ["gn.a", "gn.b"]
    gs = _restrict(gs, gns, gns)
    cof_evidence = {(a, b) for a, b in zip(gs["gn.a"], gs["gn.b"]) if a in tg_set and b in tg_set}

    net = read_table(d / "AraNet.v2.txt")
    net.columns = ["gn.a", "gn.b", "value"]
    net = _restrict(net, tg_set, tg_set)
    if net.duplicated(["gn.a", "gn.b"]).any():
        raise ValueError("duplicate AraNet pairs: acast would count them")
    net = pd.DataFrame({"row": net["gn.a"], "col": net["gn.b"], "value": net["value"].astype(float)})
    net = _rescale(net[net["value"] != 0], len(tgs) ** 2)

    # validation: SUBA3 experimental subcellular co-localization (8 compartments)
    suba = pd.read_csv(d / "Suba3-2015-2-10_7-18.csv", dtype=str, keep_default_na=False, encoding="latin-1")
    suba.columns = [re.sub(r"[^A-Za-z0-9._]", ".", c) for c in suba.columns]   # R make.names
    suba["AGI.identifier"] = suba["AGI.identifier"].str.replace(r"\..*", "", regex=True)
    keep_comp = {"mitochondrion", "nucleus", "plastid", "peroxisome", "extracellular", "endoplasmic reticulum", "vacuole", "golgi"}
    g2c: dict[str, set[str]] = {}
    for _, r in suba[suba["AGI.identifier"].isin(tg_set)].iterrows():
        comps = set()
        for c in ("location_gfp", "location_swissprot", "location_ms.ms", "location_tair"):
            comps |= {s.split(":")[0] for s in r[c].split(";") if s != ""}
        g2c.setdefault(r["AGI.identifier"], set()).update(comps & keep_comp)
    suba_pairs = _colocalization(g2c, tg_set)

    validation = {"AtRegNet": atregnet, "SUBA_colocalization": suba_pairs}
    aracyc = d / "aracyc-12.0.1-pathways.20150407"
    if aracyc.exists():
        ac = read_table(aracyc, header=True)
        ac = ac[ac["Gene-id"] != "unknown"] if "Gene-id" in ac else ac[ac["Gene.id"] != "unknown"]
        gcol = "Gene-id" if "Gene-id" in ac else "Gene.id"
        pcol = "Pathway-id" if "Pathway-id" in ac else "Pathway.id"
        pairs: Pairs = set()
        for _, genes in ac.groupby(pcol)[gcol]:
            g = sorted(set(genes) & tg_set)
            if len(g) >= 2:
                pairs |= {(a, b) for a in g for b in g if a != b}
        validation["AraCyc_pathways"] = pairs
    else:
        warnings.warn(f"AraCyc pathway file not found - skipping metabolic pathway validation: {aracyc}")

    rate = read_rds_dataframe(d / "df.rate_density_de_jaccard_coreg_tfs_1_1000.rds")
    symbols, families = {}, {}
    sym = d / "Locus_Primary_Gene_Symbol_20130117.txt"
    if sym.exists():
        s = read_table(sym, header=True)
        s = s[s["primary_gene_symbol"].fillna("") != ""].drop_duplicates("locus_name")
        symbols = dict(zip(s["locus_name"].str.upper(), s["primary_gene_symbol"]))
    fam = d / "Ath_tf_ptfb.txt"
    if fam.exists():
        f = read_table(fam, header=True)
        families = dict(zip(f["gene_model"].str.replace(r"\..*", "", regex=True).str.upper(), f["family"]))
    return Dataset("athaliana", tfs, tgs, grn, net, reg_evidence, cof_evidence, rate, validation, symbols, families)


# ---------------------------------------------------------------- D. melanogaster
def load_drosophila(root: Path) -> Dataset:
    d = Path(root) / "drosophila_datasets"
    redfly = read_table(d / "redfly.txt")
    redfly.columns = ["TF", "Target", "V3"]
    redfly = redfly[pd.to_numeric(redfly["V3"]) != 0]

    values, tfs, tgs, sparse = read_rds_matrix(d / "mat.grn.rapid.drosophila.rds")
    grn = matrix_long(values, tfs, tgs, sparse)
    gns = set(tfs) | set(tgs)
    tf_set, tg_set = set(tfs), set(tgs)
    reg_evidence = {(a, b) for a, b in acast_pairs(redfly, "TF", "Target") if a in tf_set and b in tg_set}

    gs = read_table(d / "FlyNet.GS.txt")
    gs.columns = ["gn.a", "gn.b"]
    gs = _restrict(gs, gns, gns)
    cof_evidence = {(a, b) for a, b in zip(gs["gn.a"], gs["gn.b"]) if a in tg_set and b in tg_set}

    net = read_table(d / "FlyNet.txt")
    net.columns = ["gn.a", "gn.b", "value"]
    net = _restrict(net, tg_set, tg_set)
    if net.duplicated(["gn.a", "gn.b"]).any():
        raise ValueError("duplicate FlyNet pairs: acast would count them")
    net = pd.DataFrame({"row": net["gn.a"], "col": net["gn.b"], "value": net["value"].astype(float)})
    net = _rescale(net[net["value"] != 0], len(tgs) ** 2)

    # validation sets (presence after acast + != 0 + rescaling)
    hic = read_table(d / "long_range_interactions.txt")
    hic.columns = ["gn.a", "gn.b", "value"]
    hic = _restrict(hic, tg_set, tg_set)
    chip = read_table(d / "chip_net.txt")
    chip.columns = ["gn.a", "gn.b", "value"]
    chip = _restrict(chip, tf_set, tg_set)
    imago = read_table(d / "imago.txt")
    g2c: dict[str, set[str]] = {}
    for g, t in zip(imago["V1"], imago["V2"]):
        g2c.setdefault(g, set()).add(t)
    validation = {"HiC": acast_pairs(hic, "gn.a", "gn.b"), "ChIP": acast_pairs(chip, "gn.a", "gn.b"),
                  "ImaGO_colocalization": _colocalization(g2c, tg_set)}

    rate = read_rds_dataframe(d / "df.rate_density_GENIE3_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    return Dataset("drosophila", tfs, tgs, grn, net, reg_evidence, cof_evidence, rate, validation)


def read_rds_dataframe(path: Path) -> pd.DataFrame:
    import pyreadr

    return pyreadr.read_r(path)[None]


LOADERS = {"athaliana": load_athaliana, "drosophila": load_drosophila}


def load_dataset(species: str, root: Path) -> Dataset:
    """root: folder containing the unpacked athaliana_datasets/ and drosophila_datasets/ folders."""
    return LOADERS[species](Path(root))
