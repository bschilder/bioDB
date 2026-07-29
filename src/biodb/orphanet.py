"""Orphanet rare-disease epidemiology (prevalence) client.

Wraps Orphadata's free, **CC-BY-4.0** ``en_product9_prev`` product — rare-disease
prevalence per ORPHAcode — served as a single XML file at
``www.orphadata.com``. The immediate consumer is liability-threshold effect
modeling for no-GWAS rare diseases (seq2gwas): a per-disease prevalence ``K``
sets the liability threshold ``T = Φ⁻¹(1 − K)``.

Orphanet reports prevalence as a **class band** (e.g. ``"1-9 / 100 000"``) per
``(disorder, PrevalenceType, geography, validation)`` estimate. We keep the
estimates usable as a point prevalence (``Point prevalence`` ≫ ``Prevalence at
birth`` ≫ ``Lifetime Prevalence``; **not** ``Cases/families`` counts or
``Annual incidence`` rates), prefer ``Validated`` + ``Worldwide``, and map the
class to a per-individual K midpoint (:data:`PREVALENCE_CLASS_K`).

Examples
--------
>>> from biodb.orphanet import load_prevalence, prevalence_map   # doctest: +SKIP
>>> prev = load_prevalence()                                     # doctest: +SKIP
>>> k = prevalence_map()                                         # doctest: +SKIP
>>> k["558"]  # Marfan syndrome ORPHAcode → K                    # doctest: +SKIP
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from xml.etree import ElementTree as ET

import polars as pl
import requests

logger = logging.getLogger(__name__)

PRODUCT_URL = "https://www.orphadata.com/data/xml/en_product9_prev.xml"
"""Orphadata ``en_product9_prev`` (rare-disease epidemiology / prevalence), CC-BY-4.0."""

CACHE_DIR = Path("~/.cache/biodb/orphanet").expanduser()
CACHE_DIR.mkdir(parents=True, exist_ok=True)

PREVALENCE_CLASS_K: dict[str, float] = {
    ">1 / 1000": 2e-3,
    "6-9 / 10 000": 7.5e-4,
    "1-5 / 10 000": 3e-4,
    "1-9 / 100 000": 5e-5,
    "1-9 / 1 000 000": 5e-6,
    "<1 / 1 000 000": 5e-7,
}
"""Orphanet prevalence-class band → per-individual K (band midpoint).

``"Unknown"`` / ``"Not yet documented"`` / empty are intentionally absent → no K.
"""

# Lower rank = preferred as a point-prevalence estimate. Types absent here
# (``Cases/families``, ``Annual incidence``) are excluded — they are not a
# per-individual prevalence.
_TYPE_RANK: dict[str, int] = {
    "Point prevalence": 0,
    "Prevalence at birth": 1,
    "Lifetime Prevalence": 2,
}

_USER_AGENT = "biodb/0.1 (+https://github.com/bschilder/bioDB)"
_DEFAULT_TIMEOUT = 120
_RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})


def download_prevalence(
    *, force: bool = False, max_retries: int = 5, session: requests.Session | None = None
) -> Path:
    """Download the ``en_product9_prev`` XML to the cache and return its path.

    Retries transient HTTP statuses (429, 5xx) with exponential backoff.

    Parameters
    ----------
    force : bool, default False
        Re-download even if cached.
    max_retries : int, default 5
        Maximum attempts on transient failures.
    session : requests.Session, optional
        Reuse an existing session.

    Returns
    -------
    pathlib.Path
        Path to the cached XML.
    """
    dst = CACHE_DIR / "en_product9_prev.xml"
    if dst.exists() and not force:
        return dst
    sess = session or requests.Session()
    headers = {"User-Agent": _USER_AGENT}
    delay = 1.0
    last_exc: Exception | None = None
    for attempt in range(max_retries):
        try:
            resp = sess.get(PRODUCT_URL, headers=headers, timeout=_DEFAULT_TIMEOUT)
            if resp.status_code in _RETRY_STATUSES:
                time.sleep(delay)
                delay = min(delay * 2, 30)
                continue
            resp.raise_for_status()
            dst.write_bytes(resp.content)
            logger.info("Downloaded Orphanet prevalence → %s (%d bytes)", dst, len(resp.content))
            return dst
        except requests.RequestException as exc:
            last_exc = exc
            logger.debug("Orphanet download error (attempt %d): %s", attempt + 1, exc)
            time.sleep(delay)
            delay = min(delay * 2, 30)
    raise RuntimeError(
        f"Failed to download {PRODUCT_URL} after {max_retries} attempts"
    ) from last_exc


def _en_name(parent: ET.Element, tag: str) -> str | None:
    """Return the English ``<Name>`` text of ``parent``'s ``tag`` child, if any."""
    child = parent.find(tag)
    if child is None:
        return None
    name = child.find("Name")
    return name.text.strip() if name is not None and name.text else None


def parse_prevalence(path: str | Path | None = None) -> pl.DataFrame:
    """Parse the ``en_product9_prev`` XML into a long prevalence table.

    One row per ``<Prevalence>`` estimate that (a) is a usable prevalence type
    (:data:`_TYPE_RANK`) and (b) has a class mappable via
    :data:`PREVALENCE_CLASS_K`. Streams the document (``iterparse``) so the 16 MB
    file never fully materializes as a tree.

    Parameters
    ----------
    path : str or pathlib.Path, optional
        XML path. Defaults to the cached download (call
        :func:`download_prevalence` first, or use :func:`load_prevalence`).

    Returns
    -------
    polars.DataFrame
        Columns: ``orpha_code, name, prevalence_type, prevalence_class,
        geographic, validation_status, val_moy, k``.
    """
    xml_path = Path(path) if path is not None else CACHE_DIR / "en_product9_prev.xml"
    rows: list[dict] = []
    for _event, disorder in ET.iterparse(str(xml_path), events=("end",)):
        if disorder.tag != "Disorder":
            continue
        code_el = disorder.find("OrphaCode")
        if code_el is None or not code_el.text:
            disorder.clear()
            continue
        orpha_code = code_el.text.strip()
        name_el = disorder.find("Name")
        name = name_el.text.strip() if name_el is not None and name_el.text else None
        plist = disorder.find("PrevalenceList")
        if plist is not None:
            for prev in plist.findall("Prevalence"):
                ptype = _en_name(prev, "PrevalenceType")
                if ptype not in _TYPE_RANK:
                    continue
                pclass = _en_name(prev, "PrevalenceClass")
                k = PREVALENCE_CLASS_K.get(pclass or "")
                if k is None:
                    continue
                valmoy_el = prev.find("ValMoy")
                try:
                    val_moy = (
                        float(valmoy_el.text) if valmoy_el is not None and valmoy_el.text else None
                    )
                except ValueError:
                    val_moy = None
                rows.append(
                    {
                        "orpha_code": orpha_code,
                        "name": name,
                        "prevalence_type": ptype,
                        "prevalence_class": pclass,
                        "geographic": _en_name(prev, "PrevalenceGeographic"),
                        "validation_status": _en_name(prev, "PrevalenceValidationStatus"),
                        "val_moy": val_moy,
                        "k": k,
                    }
                )
        disorder.clear()
    return pl.DataFrame(
        rows,
        schema={
            "orpha_code": pl.Utf8,
            "name": pl.Utf8,
            "prevalence_type": pl.Utf8,
            "prevalence_class": pl.Utf8,
            "geographic": pl.Utf8,
            "validation_status": pl.Utf8,
            "val_moy": pl.Float64,
            "k": pl.Float64,
        },
    )


def load_prevalence(*, force: bool = False) -> pl.DataFrame:
    """Download (if needed), parse, and cache the prevalence table as parquet.

    Parameters
    ----------
    force : bool, default False
        Re-download + re-parse even if the parquet cache exists.

    Returns
    -------
    polars.DataFrame
        The long prevalence table (see :func:`parse_prevalence`).
    """
    cache = CACHE_DIR / "prevalence.parquet"
    if cache.exists() and not force:
        return pl.read_parquet(cache)
    xml_path = download_prevalence(force=force)
    df = parse_prevalence(xml_path)
    df.write_parquet(cache)
    return df


def prevalence_map(
    prevalence: pl.DataFrame | None = None,
    *,
    prefer_worldwide: bool = True,
    validated_only: bool = False,
) -> dict[str, float]:
    """Reduce the long prevalence table to one best ``K`` per ORPHAcode.

    Selection per disorder: usable prevalence type of lowest :data:`_TYPE_RANK`,
    then ``Validated`` over ``Not yet validated``, then ``Worldwide`` geography
    (when ``prefer_worldwide``). Ties broken by the rarer (smaller) ``K`` — the
    conservative choice for a threshold model.

    Parameters
    ----------
    prevalence : polars.DataFrame, optional
        A parsed table (see :func:`parse_prevalence`); loaded via
        :func:`load_prevalence` when omitted.
    prefer_worldwide : bool, default True
        Rank ``Worldwide`` estimates above region-specific ones.
    validated_only : bool, default False
        Drop ``Not yet validated`` estimates entirely.

    Returns
    -------
    dict of str to float
        ``ORPHAcode → K`` (best estimate).
    """
    df = prevalence if prevalence is not None else load_prevalence()
    if validated_only:
        df = df.filter(pl.col("validation_status") == "Validated")
    if df.is_empty():
        return {}
    grank = (
        (pl.col("geographic") != "Worldwide").cast(pl.Int8)
        if prefer_worldwide
        else pl.lit(0, dtype=pl.Int8)
    )
    ranked = df.with_columns(
        pl.col("prevalence_type")
        .replace_strict(_TYPE_RANK, default=99, return_dtype=pl.Int64)
        .alias("_trank"),
        (pl.col("validation_status") != "Validated").cast(pl.Int8).alias("_vrank"),
        grank.alias("_grank"),
    ).with_columns(
        (pl.col("_trank") * 100 + pl.col("_vrank") * 10 + pl.col("_grank")).alias("_rank")
    )
    # Keep each disorder's best-preference tier, then take the MEDIAN K within it —
    # neutral across conflicting bands (a "smallest-K" tie-break systematically
    # over-states rarity and inflates the downstream liability effect).
    best_rank = ranked.group_by("orpha_code").agg(pl.col("_rank").min().alias("_best_rank"))
    best = (
        ranked.join(best_rank, on="orpha_code")
        .filter(pl.col("_rank") == pl.col("_best_rank"))
        .group_by("orpha_code")
        .agg(pl.col("k").median().alias("k"))
    )
    return dict(best.select("orpha_code", "k").iter_rows())
