"""GRACE command line.

    python -m grace reproduce --data DIR [species ...]     published models -> published networks
    python -m grace train     --data DIR species [--models 100] [--workers N] [--seed-offset K] [--out DIR]
    python -m grace validate  --data DIR species [--network published | path/to/network.csv]

DIR holds the unpacked tutorial datasets (athaliana_datasets/, drosophila_datasets/); it can also be set
with the GRACE_DATA_DIR environment variable.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from .pipeline import ensemble_network, prepare
from .reference import SPECIES, SPECIES_NAMES, published_models, published_network
from .train import default_workers, train_models
from .validation import evidence_sets, validate


def _data(args) -> Path:
    d = args.data or os.environ.get("GRACE_DATA_DIR")
    if not d:
        sys.exit("error: pass --data DIR or set GRACE_DATA_DIR to the folder holding the tutorial datasets")
    return Path(d)


def _overlap(a: set, b: set) -> str:
    j = len(a & b) / len(a | b) if a | b else 1.0
    return f"{len(a & b)} shared, Jaccard {j:.3f}"


def cmd_reproduce(args):
    ok = True
    for s in args.species or SPECIES:
        t = time.perf_counter()
        p = prepare(s, _data(args))
        net = ensemble_network(p, published_models(s))
        got, ref = set(zip(net.TF, net.TG)), published_network(s)
        same = got == ref
        ok &= same
        print(f"{SPECIES_NAMES[s]}: {len(got)} links vs {len(ref)} published -> "
              f"{'identical' if same else 'DIFFERENT (' + _overlap(got, ref) + ')'}  [{time.perf_counter() - t:.1f} s]")
    return 0 if ok else 1


def cmd_train(args):
    p = prepare(args.species, _data(args))
    seeds = [1234 + 158 * i + args.seed_offset for i in range(1, args.models + 1)]
    t = time.perf_counter()
    models = train_models(p.training, n_workers=args.workers, seeds=seeds)
    dt = time.perf_counter() - t
    net = ensemble_network(p, models, threshold=args.threshold)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = f"{args.species}" + (f"_seed{args.seed_offset}" if args.seed_offset else "")
    pd.DataFrame([asdict(m) for m in models]).to_csv(out / f"models_{tag}.csv", index=False)
    net.to_csv(out / f"network_{tag}.csv", index=False)
    print(f"{SPECIES_NAMES[args.species]}: {len(models)} models in {dt:.1f} s on {args.workers} workers")
    print(f"ensemble network: {len(net)} links (support > {args.threshold}) -> {out / f'network_{tag}.csv'}")
    print(f"vs published network: {_overlap(set(zip(net.TF, net.TG)), published_network(args.species))}")
    return 0


def cmd_validate(args):
    p = prepare(args.species, _data(args))
    if args.network == "published":
        net = published_network(args.species)
    else:
        df = pd.read_csv(args.network)
        net = set(zip(df["TF"], df["TG"]))
    links, pairs = evidence_sets(p.dataset)
    res = validate(net, p.network.tfs, p.network.tgs, links, pairs)
    print(f"{SPECIES_NAMES[args.species]}: {len(net)} links")
    print(f"{'evidence':<40}{'kind':<12}{'predicted':>10}{'hits':>8}{'fold':>8}{'p':>11}")
    for r in res:
        print(f"{r.evidence:<40}{r.kind:<12}{r.n_predicted:>10}{r.hits:>8}{r.fold_change:>8.2f}{r.p_value:>11.2e}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m grace", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--data", help="folder with the unpacked tutorial datasets (or GRACE_DATA_DIR)")

    r = sub.add_parser("reproduce", parents=[common], help="rebuild the published networks from the published models")
    r.add_argument("species", nargs="*", choices=SPECIES)
    r.set_defaults(func=cmd_reproduce)

    t = sub.add_parser("train", parents=[common], help="train models and build an ensemble network")
    t.add_argument("species", choices=SPECIES)
    t.add_argument("--models", type=int, default=100)
    t.add_argument("--workers", type=int, default=default_workers(), help="parallel workers (default: cores - 2)")
    t.add_argument("--seed-offset", type=int, default=0, help="0 = the seeds of the R tutorials")
    t.add_argument("--threshold", type=float, default=0.5, help="minimum fraction of models predicting a link")
    t.add_argument("--out", default="results")
    t.set_defaults(func=cmd_train)

    v = sub.add_parser("validate", parents=[common], help="enrichment of a network for independent evidence")
    v.add_argument("species", choices=SPECIES)
    v.add_argument("--network", default="published", help="'published' or a CSV with TF and TG columns")
    v.set_defaults(func=cmd_validate)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
