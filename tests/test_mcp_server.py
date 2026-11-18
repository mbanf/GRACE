"""GRACE MCP server tools (need GRACE_DATA_DIR; the tool list check does not)."""

import asyncio
import os
from pathlib import Path

import pytest

pytest.importorskip("mcp")
from grace_mcp import server as s  # noqa: E402

DATA = Path(os.environ.get("GRACE_DATA_DIR", "data"))
HAS_DATA = (DATA / "athaliana_datasets").exists() and (DATA / "drosophila_datasets").exists()
needs_data = pytest.mark.skipif(not HAS_DATA, reason="set GRACE_DATA_DIR to the tutorial datasets")

TOOLS = {"grace_summary", "grace_query_tf", "grace_query_gene", "grace_query_link", "grace_build_network",
         "grace_train", "grace_validate", "grace_compare"}


def test_tools_registered():
    tools = asyncio.run(s.mcp.list_tools())
    assert {t.name for t in tools} == TOOLS


@pytest.fixture(scope="module")
def state():
    s.STATE.data_dir = DATA
    return s.STATE


@needs_data
def test_summary(state):
    r = s.grace_summary()
    assert r["athaliana"]["candidate_links"] == 10098
    assert r["athaliana"]["networks"]["published"] == 792
    assert r["drosophila"]["networks"]["published"] == 7164


@needs_data
def test_query_tf_by_symbol(state):
    r = s.grace_query_tf("FLC")
    assert r["tf"]["id"] == "AT5G10140" and r["tf"]["tf_family"] == "MIKC"
    assert r["targets_in_network"] == len(r["targets"]) >= 1


@needs_data
def test_query_gene_and_link(state):
    r = s.grace_query_gene("AT4G10020")
    assert any(x["id"] == "AT5G10140" for x in r["regulators"])
    link = s.grace_query_link("FLC", "AT4G10020")
    assert link["candidate"] and link["support_by_network"]["published"] > 0.5


@needs_data
def test_build_network_threshold(state):
    r = s.grace_build_network("athaliana", threshold=0.7)
    assert r["vs_published"]["only_a"] == 0 and r["links"] < 792    # stricter threshold: a subset


@needs_data
def test_validate(state):
    r = s.grace_validate("drosophila")
    chip = next(e for e in r["enrichment"] if e["evidence"] == "ChIP")
    assert chip["hits"] == 736 and chip["fold_change"] > 2
