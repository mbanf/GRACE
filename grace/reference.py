"""Published GRACE v1.0 results shipped with the package: the 100 model weights per species
(alpha, gamma, lambda) and the published ensemble networks."""

from __future__ import annotations

from importlib.resources import files

import pandas as pd

from .train import Model

SPECIES = ("athaliana", "drosophila")
SPECIES_NAMES = {"athaliana": "Arabidopsis thaliana", "drosophila": "Drosophila melanogaster"}


def _path(name: str):
    return files("grace") / "reference" / name


def published_models(species: str) -> list[Model]:
    """The 100 models of the paper (only their weights are stored)."""
    w = pd.read_csv(_path(f"published_models_{species}.tsv"), sep="\t")
    return [Model(a, g, l, fscore=float("nan"), n_pred=float("nan"), n_gs=float("nan"),
                  n_coreg_bp=float("nan"), seed=-1) for a, g, l in w[["alpha", "gamma", "lambda"]].to_numpy()]


def published_network(species: str) -> set[tuple[str, str]]:
    df = pd.read_csv(_path(f"published_network_{species}.csv"))
    return set(zip(df["TF"], df["TG"]))
