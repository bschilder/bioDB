"""Tests for :mod:`biodb.orphanet` (Orphadata ``en_product9_prev`` prevalence).

Offline tests mock the HTTP download (:mod:`responses`) and parse a small XML
fixture. The live endpoint is exercised separately (``@pytest.mark.network``).
"""

from __future__ import annotations

import pytest
import responses

from biodb import orphanet

_FIXTURE_XML = """<?xml version="1.0" encoding="UTF-8"?>
<JDBOR>
  <DisorderList count="2">
    <Disorder id="1">
      <OrphaCode>558</OrphaCode>
      <Name lang="en">Marfan syndrome</Name>
      <PrevalenceList count="2">
        <Prevalence id="1">
          <PrevalenceType id="a"><Name lang="en">Point prevalence</Name></PrevalenceType>
          <PrevalenceClass id="b"><Name lang="en">1-5 / 10 000</Name></PrevalenceClass>
          <PrevalenceGeographic id="c"><Name lang="en">Worldwide</Name></PrevalenceGeographic>
          <PrevalenceValidationStatus id="d"><Name lang="en">Validated</Name></PrevalenceValidationStatus>
          <ValMoy>3.0</ValMoy>
        </Prevalence>
        <Prevalence id="2">
          <PrevalenceType id="a"><Name lang="en">Cases/families</Name></PrevalenceType>
          <PrevalenceClass/>
          <ValMoy>10.0</ValMoy>
        </Prevalence>
      </PrevalenceList>
    </Disorder>
    <Disorder id="2">
      <OrphaCode>999</OrphaCode>
      <Name lang="en">Test disease</Name>
      <PrevalenceList count="2">
        <Prevalence id="3">
          <PrevalenceType id="a"><Name lang="en">Point prevalence</Name></PrevalenceType>
          <PrevalenceClass><Name lang="en">1-9 / 1 000 000</Name></PrevalenceClass>
          <PrevalenceGeographic><Name lang="en">Worldwide</Name></PrevalenceGeographic>
          <PrevalenceValidationStatus><Name lang="en">Validated</Name></PrevalenceValidationStatus>
        </Prevalence>
        <Prevalence id="4">
          <PrevalenceType id="a"><Name lang="en">Point prevalence</Name></PrevalenceType>
          <PrevalenceClass><Name lang="en">1-9 / 100 000</Name></PrevalenceClass>
          <PrevalenceGeographic><Name lang="en">Worldwide</Name></PrevalenceGeographic>
          <PrevalenceValidationStatus><Name lang="en">Validated</Name></PrevalenceValidationStatus>
        </Prevalence>
      </PrevalenceList>
    </Disorder>
  </DisorderList>
</JDBOR>
"""


def _write_fixture(tmp_path):  # noqa: ANN001, ANN202
    p = tmp_path / "en_product9_prev.xml"
    p.write_text(_FIXTURE_XML)
    return p


def test_class_k_map_bands() -> None:
    # The band midpoints are monotone in rarity.
    k = orphanet.PREVALENCE_CLASS_K
    assert k["1-5 / 10 000"] > k["1-9 / 100 000"] > k["1-9 / 1 000 000"] > k["<1 / 1 000 000"]
    assert "Unknown" not in k and "Not yet documented" not in k


def test_parse_prevalence_filters_usable_types_and_classes(tmp_path) -> None:  # noqa: ANN001
    df = orphanet.parse_prevalence(_write_fixture(tmp_path))
    # Marfan: only the Point-prevalence row survives (Cases/families dropped).
    marfan = df.filter(pl_col_eq(df, "orpha_code", "558"))
    assert marfan.height == 1
    assert marfan["prevalence_class"][0] == "1-5 / 10 000"
    assert marfan["k"][0] == pytest.approx(3e-4)
    # Disorder 999: both point-prevalence rows are kept.
    assert df.filter(pl_col_eq(df, "orpha_code", "999")).height == 2


def test_prevalence_map_median_tiebreak(tmp_path) -> None:  # noqa: ANN001
    df = orphanet.parse_prevalence(_write_fixture(tmp_path))
    kmap = orphanet.prevalence_map(df)
    assert kmap["558"] == pytest.approx(3e-4)
    # 999 has two worldwide+validated+point estimates (5e-6, 5e-5) → median.
    assert kmap["999"] == pytest.approx((5e-6 + 5e-5) / 2)


@responses.activate
def test_download_prevalence_caches(tmp_path, monkeypatch) -> None:  # noqa: ANN001
    monkeypatch.setattr(orphanet, "CACHE_DIR", tmp_path)
    responses.add(responses.GET, orphanet.PRODUCT_URL, body=_FIXTURE_XML, status=200)
    p = orphanet.download_prevalence()
    assert p.exists() and p.read_text().startswith("<?xml")
    # Second call is a cache hit — no further HTTP request registered/needed.
    p2 = orphanet.download_prevalence()
    assert p2 == p
    assert len(responses.calls) == 1


@responses.activate
def test_download_prevalence_retries_5xx(tmp_path, monkeypatch) -> None:  # noqa: ANN001
    monkeypatch.setattr(orphanet, "CACHE_DIR", tmp_path)
    responses.add(responses.GET, orphanet.PRODUCT_URL, status=503)
    responses.add(responses.GET, orphanet.PRODUCT_URL, body=_FIXTURE_XML, status=200)
    p = orphanet.download_prevalence(max_retries=3)
    assert p.exists()
    assert len(responses.calls) == 2


def pl_col_eq(df, col, val):  # noqa: ANN001, ANN201
    """Small helper: boolean mask ``df[col] == val`` (import-light)."""
    import polars as pl

    return pl.col(col) == val
