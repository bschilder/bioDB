"""Tests for :mod:`biodb.panukbb` — Pan-UK Biobank full GWAS summary statistics.

The pure functions (bgz decompress, manifest parse, trait selection, URL
derivation, sumstats-line parse) are exercised against **real-data fixtures**
sliced from the live Pan-UKBB release, so the offline suite is deterministic yet
faithful to the upstream schema. ``@pytest.mark.network`` tests at the bottom
probe the live S3 bucket + htslib remote tabix and are skipped in CI.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pytest

from biodb import panukbb
from tests.conftest import is_upstream_outage

FIXTURES = Path(__file__).parent / "fixtures"
_MANIFEST_TSV = FIXTURES / "panukbb_manifest_slice.tsv"
_MANIFEST_BGZ = FIXTURES / "panukbb_manifest_slice.tsv.gz"  # bgzip'd copy (valid bgzf)
_SUMSTATS_TSV = FIXTURES / "panukbb_sumstats_apoe.tsv"  # header + real APOE region lines


# --------------------------------------------------------------------------- #
# _decompress_bgz                                                             #
# --------------------------------------------------------------------------- #
def test_decompress_bgz_roundtrips_multiblock() -> None:
    """A bgzip'd (multi-member gzip) payload decompresses to the original TSV."""
    raw = _MANIFEST_BGZ.read_bytes()
    out = panukbb._decompress_bgz(raw)
    assert out == _MANIFEST_TSV.read_bytes()


def test_decompress_bgz_passes_through_plain_text() -> None:
    """Already-decompressed bytes (no gzip magic) are returned unchanged."""
    plain = b"chr\tpos\n1\t100\n"
    assert panukbb._decompress_bgz(plain) == plain


# --------------------------------------------------------------------------- #
# _parse_manifest                                                             #
# --------------------------------------------------------------------------- #
def test_parse_manifest_shape_and_columns() -> None:
    """Robust parse recovers every row/column despite embedded quotes/newlines."""
    df = panukbb._parse_manifest(_MANIFEST_TSV.read_bytes())
    assert len(df) == 15  # 15 data rows in the slice
    assert df.shape[1] == 82  # full manifest schema
    for col in ("trait_type", "phenocode", "sldsc_25bin_h2_z_EUR", "aws_path"):
        assert col in df.columns


def test_parse_manifest_accepts_bgz_bytes() -> None:
    """Passing the compressed bytes directly is transparently decompressed."""
    df = panukbb._parse_manifest(_MANIFEST_BGZ.read_bytes())
    assert len(df) == 15


def test_parse_manifest_preserves_embedded_special_chars() -> None:
    """A description field containing a quote/newline is not truncated or split.

    This is the exact failure mode that broke a naive polars ``read_csv`` — a
    quote-aware parser must keep the row intact.
    """
    df = panukbb._parse_manifest(_MANIFEST_TSV.read_bytes())
    # No row should have leaked into a phantom extra row: phenocode is always set
    assert df["phenocode"].notna().all()
    assert (df["phenocode"].astype(str).str.len() > 0).all()


# --------------------------------------------------------------------------- #
# select_traits                                                               #
# --------------------------------------------------------------------------- #
def test_select_traits_default_keeps_qc_pass_heritable_biomarkers() -> None:
    """Default filter = continuous/biomarker, EUR QC PASS, h2_z > 4.

    The slice has exactly 6 biomarker rows with QC==PASS and h2_z in [8.9, 16.6];
    every other row is either the wrong trait_type or QC != PASS.
    """
    m = panukbb._parse_manifest(_MANIFEST_TSV.read_bytes())
    sel = panukbb.select_traits(m, n=100)
    assert len(sel) == 6
    assert set(sel["trait_type"]) <= {"continuous", "biomarkers"}
    assert (sel["phenotype_qc_EUR"] == "PASS").all()
    assert (sel["sldsc_25bin_h2_z_EUR"].astype(float) > 4).all()


def test_select_traits_sorted_by_heritability_and_capped() -> None:
    """Results are ranked by h2_z (desc) and truncated to n."""
    m = panukbb._parse_manifest(_MANIFEST_TSV.read_bytes())
    sel = panukbb.select_traits(m, n=3)
    assert len(sel) == 3
    h2z = sel["sldsc_25bin_h2_z_EUR"].astype(float).tolist()
    assert h2z == sorted(h2z, reverse=True)


def test_select_traits_lambda_gc_window_is_off_by_default() -> None:
    """Polygenic high-N biomarkers have inflated lambda_gc (1.3-1.7) yet PASS QC.

    A naive lambda_gc window of [0.9, 1.2] would reject exactly the well-powered
    heritable traits — so it must be OPT-IN, not the default.
    """
    m = panukbb._parse_manifest(_MANIFEST_TSV.read_bytes())
    assert len(panukbb.select_traits(m)) == 6  # default: lambda filter off
    tight = panukbb.select_traits(m, lambda_gc_range=(0.9, 1.2))
    assert len(tight) == 0  # every PASS biomarker here is legitimately inflated


def test_select_traits_max_per_category_diversifies() -> None:
    """``max_per_category`` caps near-duplicate traits so the set spans biology.

    All 6 fixture PASS-biomarkers share one category; a cap of 2 must keep only
    the 2 highest-h2_z of them (the top-by-h2_z ranking is preserved within a
    category).
    """
    m = panukbb._parse_manifest(_MANIFEST_TSV.read_bytes())
    sel = panukbb.select_traits(m, max_per_category=2)
    assert len(sel) == 2
    assert sel["category"].nunique() == 1
    h2z = sel["sldsc_25bin_h2_z_EUR"].astype(float).tolist()
    assert h2z == sorted(h2z, reverse=True)
    # they are the two most-heritable of the six
    all6 = panukbb.select_traits(m)["sldsc_25bin_h2_z_EUR"].astype(float).nlargest(2).tolist()
    assert h2z == all6


def test_select_traits_fetches_manifest_when_none(monkeypatch) -> None:
    """Passing ``manifest=None`` pulls it via ``phenotype_manifest`` (no network)."""
    m = panukbb._parse_manifest(_MANIFEST_TSV.read_bytes())
    monkeypatch.setattr(panukbb, "phenotype_manifest", lambda force=False: m)
    sel = panukbb.select_traits(n=2)
    assert len(sel) == 2


# --------------------------------------------------------------------------- #
# _sumstats_urls                                                              #
# --------------------------------------------------------------------------- #
def test_sumstats_urls_from_aws_path_string() -> None:
    """s3:// flat-file path → https data URL + the sibling tabix-index URL."""
    aws = "s3://pan-ukb-us-east-1/sumstats_flat_files/biomarkers-30600-both_sexes-irnt.tsv.bgz"
    data_url, idx_url = panukbb._sumstats_urls(aws)
    assert data_url == (
        "https://pan-ukb-us-east-1.s3.amazonaws.com/"
        "sumstats_flat_files/biomarkers-30600-both_sexes-irnt.tsv.bgz"
    )
    assert idx_url == (
        "https://pan-ukb-us-east-1.s3.amazonaws.com/"
        "sumstats_flat_files_tabix/biomarkers-30600-both_sexes-irnt.tsv.bgz.tbi"
    )


def test_regions_to_bed_0based_halfopen_and_prefix_stripped() -> None:
    """1-based-inclusive regions → BED (0-based half-open), ``chr`` stripped."""
    bed = panukbb._regions_to_bed([("chr19", 100, 200), ("1", 5, 5)])
    assert bed == "19\t99\t200\n1\t4\t5\n"


def test_regions_to_bed_empty() -> None:
    assert panukbb._regions_to_bed([]) == ""


def test_sumstats_urls_from_manifest_row() -> None:
    """A manifest row (dict/Series) resolves via its aws_path / aws_path_tabix."""
    row = {
        "aws_path": "s3://pan-ukb-us-east-1/sumstats_flat_files/x.tsv.bgz",
        "aws_path_tabix": "s3://pan-ukb-us-east-1/sumstats_flat_files_tabix/x.tsv.bgz.tbi",
    }
    data_url, idx_url = panukbb._sumstats_urls(row)
    assert data_url.endswith("sumstats_flat_files/x.tsv.bgz")
    assert idx_url.endswith("sumstats_flat_files_tabix/x.tsv.bgz.tbi")


# --------------------------------------------------------------------------- #
# _ensure_local_index (index caching — avoids re-fetch + cwd litter)          #
# --------------------------------------------------------------------------- #
def test_ensure_local_index_caches_and_reuses(tmp_path, monkeypatch) -> None:
    """The remote ``.tbi`` is downloaded once, then served from the local cache.

    htslib otherwise re-fetches the index into cwd on every region query;
    caching it locally is both a big speedup for many-window workloads and keeps
    the caller's working directory clean.
    """
    monkeypatch.setattr(panukbb, "CACHE_DIR", tmp_path)
    calls = {"n": 0}

    class _FakeResp:
        content = b"TBI-BYTES"

        def raise_for_status(self) -> None:
            pass

    class _FakeSession:
        def get(self, url, **kw):  # noqa: ANN001
            calls["n"] += 1
            return _FakeResp()

    monkeypatch.setattr(panukbb, "_session", lambda: _FakeSession())

    url = "https://x/sumstats_flat_files_tabix/foo.tsv.bgz.tbi"
    p1 = panukbb._ensure_local_index(url)
    p2 = panukbb._ensure_local_index(url)
    assert p1 == p2
    assert p1.read_bytes() == b"TBI-BYTES"
    assert calls["n"] == 1  # second call is a cache hit — no re-download


# --------------------------------------------------------------------------- #
# _parse_sumstats_lines                                                       #
# --------------------------------------------------------------------------- #
def _sumstats_fixture_header_and_lines() -> tuple[list[str], list[str]]:
    text = _SUMSTATS_TSV.read_text().splitlines()
    return text[0].split("\t"), text[1:]


def test_parse_sumstats_lines_tidy_schema() -> None:
    """Region lines → tidy per-variant frame with chrom/pos/ref/alt + EUR stats."""
    header, lines = _sumstats_fixture_header_and_lines()
    df = panukbb._parse_sumstats_lines(header, lines, ancestries=("EUR",))
    assert list(df.columns[:4]) == ["chrom", "pos", "ref", "alt"]
    assert {"af_EUR", "beta_EUR", "se_EUR", "neglog10_pval_EUR"} <= set(df.columns)
    assert len(df) == len(lines)
    assert df["pos"].dtype.kind in "iu"  # integer positions
    assert df["beta_EUR"].dtype.kind == "f"  # numeric effect sizes


def test_parse_sumstats_lines_computes_z() -> None:
    """z = beta/se is derived per ancestry when requested."""
    header, lines = _sumstats_fixture_header_and_lines()
    df = panukbb._parse_sumstats_lines(header, lines, ancestries=("EUR",), add_z=True)
    assert "z_EUR" in df.columns
    expected = df["beta_EUR"] / df["se_EUR"]
    pd.testing.assert_series_equal(df["z_EUR"], expected, check_names=False)


def test_parse_sumstats_lines_unknown_ancestry_raises() -> None:
    """Requesting an ancestry the file lacks fails loudly, not silently."""
    header, lines = _sumstats_fixture_header_and_lines()
    with pytest.raises(KeyError, match="ZZZ"):
        panukbb._parse_sumstats_lines(header, lines, ancestries=("ZZZ",))


# --------------------------------------------------------------------------- #
# Live integration — gated behind --run-network                              #
# --------------------------------------------------------------------------- #
@pytest.mark.network
def test_live_phenotype_manifest() -> None:
    try:
        m = panukbb.phenotype_manifest()
    except Exception as exc:  # noqa: BLE001
        if is_upstream_outage(exc):
            pytest.skip(f"Pan-UKBB S3 unavailable: {exc}")
        raise
    assert len(m) > 5000
    assert "aws_path" in m.columns


@pytest.mark.network
@pytest.mark.slow
@pytest.mark.skipif(shutil.which("tabix") is None, reason="htslib 'tabix' binary not installed")
def test_live_sumstats_region_apoe() -> None:
    aws = "s3://pan-ukb-us-east-1/sumstats_flat_files/biomarkers-30600-both_sexes-irnt.tsv.bgz"
    try:
        df = panukbb.sumstats_region(aws, "19", 45405000, 45415000, ancestries=("EUR",))
    except Exception as exc:  # noqa: BLE001
        if is_upstream_outage(exc):
            pytest.skip(f"Pan-UKBB S3 unavailable: {exc}")
        raise
    assert len(df) > 50  # dense sumstats — every variant, not just genome-wide-sig
    assert df["beta_EUR"].notna().any()


@pytest.mark.network
@pytest.mark.slow
@pytest.mark.skipif(shutil.which("tabix") is None, reason="htslib 'tabix' binary not installed")
def test_live_sumstats_regions_multi() -> None:
    """One tabix -R pass returns variants from multiple windows across contigs."""
    aws = "s3://pan-ukb-us-east-1/sumstats_flat_files/biomarkers-30600-both_sexes-irnt.tsv.bgz"
    regions = [("19", 45405000, 45410000), ("1", 55505000, 55515000)]  # APOE + PCSK9 (hg19)
    try:
        df = panukbb.sumstats_regions(aws, regions, ancestries=("EUR",))
    except Exception as exc:  # noqa: BLE001
        if is_upstream_outage(exc):
            pytest.skip(f"Pan-UKBB S3 unavailable: {exc}")
        raise
    assert len(df) > 50
    assert (df["chrom"] == "19").any() and (df["chrom"] == "1").any()  # both windows present
