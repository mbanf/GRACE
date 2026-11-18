#!/usr/bin/env bash
# Run a GRACE v1.0 tutorial on current R.
#
# Usage: ./run_tutorial.sh athaliana|drosophila /path/to/datasets
#   /path/to/datasets must contain the unpacked athaliana_datasets/ and/or drosophila_datasets/ folders.
#
# The tutorials source all GRACE files from their working directory, so this script builds a flat
# working folder (work/) with the R sources and links to the datasets, then runs the tutorial there.
set -euo pipefail
species="${1:?species: athaliana or drosophila}"
data="$(cd "${2:?path to datasets}" && pwd)"
here="$(cd "$(dirname "$0")" && pwd)"
work="$here/work"
mkdir -p "$work"
cp "$here"/R/*.R "$here"/tutorials/*.R "$work"/
ln -sfn "$data/${species}_datasets" "$work/${species}_datasets"
cd "$work"
Rscript "GRACE_tutorial_${species}.R"
