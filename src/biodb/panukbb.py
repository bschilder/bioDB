"""Pan-UK Biobank (Pan-UKBB) full GWAS summary statistics.

Pan-UKBB ran a standardised multi-ancestry GWAS of ~7,200 phenotypes across
~500k UK Biobank participants and released the results as public
`AWS Open Data <https://registry.opendata.aws/broad-pan-ukb/>`_
(bucket ``s3://pan-ukb-us-east-1``, also served over https). Unlike credible-set
resources (Open Targets, hg-horizon fine-mapped), these are **full, dense
summary statistics** — every one of ~28M imputed variants gets a per-ancestry
``(af, beta, se, neglog10_pval)``, genome-wide-significant or not.

Three entry points:

- :func:`phenotype_manifest` -- the per-phenotype manifest as a tidy
  ``DataFrame`` (power, heritability, QC, per-ancestry N, flat-file paths).
- :func:`select_traits` -- pick similarly-powered, heritable, QC-passing
  phenotypes. The default (continuous + biomarker, EUR QC ``PASS``, ``h2_z``
  thresholded) is *power-matched by construction*: those traits are measured on
  ~the full EUR cohort, so their effective N is near-uniform.
- :func:`sumstats_region` -- htslib-tabix a genomic region of one phenotype's
  flat file into a tidy per-variant ``DataFrame``. Streams the region over the
  https URL via HTTP range requests — **no multi-GB download**.

Notes
-----
* **Genome build is GRCh37/hg19.** Flat-file contigs are ``1..22, X`` (no
  ``chr`` prefix). Consumers aligning to GRCh38 must liftover.
* ``beta`` is oriented to **ALT = effect allele**; ``af_{anc}`` is the effect
  (ALT) allele frequency *in that ancestry's GWAS cohort* (not a reference
  panel) — the quantity needed to weight/anneal per-variant effects.
* ``sumstats_region`` needs the htslib ``tabix`` binary on ``PATH`` (default
  backend) *or* the optional ``pysam`` package (``backend="pysam"``). Every
  seqlab/bcftools environment already ships htslib.
"""

from __future__ import annotations

import csv
import gzip
import io
import logging
import os
import shutil
import subprocess
import tempfile
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import requests

logger = logging.getLogger(__name__)

S3_BUCKET = "pan-ukb-us-east-1"
S3_BASE = f"https://{S3_BUCKET}.s3.amazonaws.com"
"""HTTPS root of the public Pan-UKBB Open Data bucket."""

MANIFEST_URL = f"{S3_BASE}/sumstats_release/phenotype_manifest.tsv.bgz"
"""Per-phenotype manifest (block-gzipped TSV, ~7,200 rows)."""

DEFAULT_RELEASE = "0.4"
"""Pan-UKBB release these paths correspond to (r0.4, the public flat files)."""

CACHE_DIR = Path("~/.cache/biodb/panukbb").expanduser() / DEFAULT_RELEASE
CACHE_DIR.mkdir(parents=True, exist_ok=True)

ANCESTRY_CODES: tuple[str, ...] = ("AFR", "AMR", "CSA", "EAS", "EUR", "MID")
"""The six Pan-UKBB genetic-ancestry groups (``meta`` = inverse-variance meta-analysis)."""

_USER_AGENT = "biodb/0.1 (+https://github.com/bschilder/bioDB)"

# Per-ancestry summary-statistic columns present in every flat file (those
# ancestries that passed QC). ``low_confidence`` is a per-variant QC flag.
_PER_ANCESTRY_STATS: tuple[str, ...] = ("af", "beta", "se", "neglog10_pval", "low_confidence")


# --------------------------------------------------------------------------- #
# Decompression + manifest parsing (pure)                                     #
# --------------------------------------------------------------------------- #
def _decompress_bgz(data: bytes) -> bytes:
    """Decompress block-gzip (bgzf) bytes, or return plain bytes unchanged.

    bgzf is a valid *multi-member* gzip stream, so :class:`gzip.GzipFile` reads
    it end-to-end. Bytes lacking the gzip magic (``1f 8b``) are already plain
    text and are returned as-is — this lets callers pass either form.

    Parameters
    ----------
    data
        Raw bytes, either bgzf/gzip-compressed or plain.

    Returns
    -------
    bytes
        Decompressed (or pass-through) bytes.
    """
    if data[:2] != b"\x1f\x8b":
        return data
    with gzip.GzipFile(fileobj=io.BytesIO(data)) as fh:
        return fh.read()


def _parse_manifest(data: bytes | str) -> pd.DataFrame:
    """Parse the Pan-UKBB phenotype manifest robustly into a string ``DataFrame``.

    The manifest's ``description``/``description_more`` fields contain embedded
    quotes and newlines inside quoted cells — which breaks a naive line-oriented
    reader (this is the exact failure that makes polars' native ``read_csv`` see
    phantom rows). Python's quote-aware :mod:`csv` reader parses it faithfully.

    Parameters
    ----------
    data
        Manifest bytes (bgzf-compressed or plain) or a decoded string.

    Returns
    -------
    pandas.DataFrame
        Every column as ``str`` (raw). Numeric coercion is deferred to
        :func:`select_traits`, which knows which columns are numeric.
    """
    text = _decompress_bgz(data).decode("utf-8") if isinstance(data, bytes) else data
    reader = csv.reader(io.StringIO(text), delimiter="\t")
    rows = list(reader)
    header, body = rows[0], rows[1:]
    return pd.DataFrame(body, columns=header, dtype="string")


# --------------------------------------------------------------------------- #
# Network: fetch the manifest (cached)                                        #
# --------------------------------------------------------------------------- #
def _session() -> requests.Session:
    sess = requests.Session()
    sess.headers.update({"User-Agent": _USER_AGENT})
    return sess


def phenotype_manifest(force: bool = False) -> pd.DataFrame:
    """Fetch (and cache) the Pan-UKBB per-phenotype manifest.

    Parameters
    ----------
    force
        Re-download even if a cached copy exists.

    Returns
    -------
    pandas.DataFrame
        ~7,200 rows × 82 columns (string dtype). Key columns:
        ``trait_type, phenocode, description, category, pops_pass_qc,
        n_cases_EUR, n_controls_EUR, sldsc_25bin_h2_z_EUR, lambda_gc_EUR,
        phenotype_qc_EUR, filename, aws_path, aws_path_tabix``.
    """
    cache = CACHE_DIR / "phenotype_manifest.parquet"
    if cache.exists() and not force:
        return pd.read_parquet(cache)
    logger.info("Fetching Pan-UKBB phenotype manifest from %s", MANIFEST_URL)
    resp = _session().get(MANIFEST_URL, timeout=120)
    resp.raise_for_status()
    df = _parse_manifest(resp.content)
    df.to_parquet(cache)
    return df


# --------------------------------------------------------------------------- #
# Trait selection (pure)                                                       #
# --------------------------------------------------------------------------- #
def _h2_z_column(ancestry: str) -> str:
    """Heritability-z column for an ancestry (EUR uses S-LDSC, others RHE-mc)."""
    return "sldsc_25bin_h2_z_EUR" if ancestry == "EUR" else f"rhemc_25bin_50rv_h2_z_{ancestry}"


def select_traits(
    manifest: pd.DataFrame | None = None,
    *,
    n: int = 100,
    trait_types: tuple[str, ...] = ("continuous", "biomarkers"),
    ancestry: str = "EUR",
    min_h2_z: float = 4.0,
    require_qc_pass: bool = True,
    lambda_gc_range: tuple[float, float] | None = None,
    n_cases_range: tuple[float, float] | None = None,
    max_per_category: int | None = None,
) -> pd.DataFrame:
    """Select similarly-powered, heritable, QC-passing phenotypes.

    The default is deliberately conservative and **power-matched by
    construction**: continuous + biomarker phenotypes are measured on ~the whole
    EUR cohort, so their effective N is near-uniform; requiring Pan-UKBB's own
    ``phenotype_qc_{anc} == "PASS"`` and a confident heritability-z
    (``sldsc_25bin_h2_z > min_h2_z``) keeps only traits with real, well-estimated
    genetic signal.

    ``lambda_gc`` filtering is **off by default on purpose**. For high-N
    polygenic biomarkers, genomic inflation ``lambda_gc`` is legitimately
    elevated (1.3-1.7) by true polygenicity, not confounding — a naive
    ``[0.9, 1.2]`` window would reject exactly the best-powered heritable traits.
    Pan-UKBB's curated ``phenotype_qc`` flag already accounts for this, so trust
    it and leave the window opt-in.

    Parameters
    ----------
    manifest
        A manifest ``DataFrame`` (from :func:`phenotype_manifest`). If ``None``,
        it is fetched.
    n
        Keep the top-``n`` by heritability-z.
    trait_types
        Which ``trait_type`` values to keep.
    ancestry
        Ancestry whose QC / heritability / N columns gate the selection.
    min_h2_z
        Minimum ``sldsc_25bin_h2_z_{ancestry}`` (confidence that h2 > 0).
    require_qc_pass
        Require ``phenotype_qc_{ancestry} == "PASS"``.
    lambda_gc_range
        Optional ``(lo, hi)`` window on ``lambda_gc_{ancestry}`` (opt-in; see
        above).
    n_cases_range
        Optional ``(lo, hi)`` window on ``n_cases_{ancestry}`` (tighten the
        power match beyond the trait-type heuristic).
    max_per_category
        Optional cap on how many traits to keep per ``category`` value, applied
        after ranking (keeps the most-heritable within each category). Promotes
        biological diversity — e.g. it prevents the four near-identical
        "leg/arm fat" impedance traits from crowding out the rest. Rows with a
        missing/``"NA"`` category are exempt (each is kept).

    Returns
    -------
    pandas.DataFrame
        The selected rows (original columns), ranked by heritability-z
        descending, truncated to ``n``.
    """
    if manifest is None:
        manifest = phenotype_manifest()
    m = manifest.copy()

    h2z_col = _h2_z_column(ancestry)
    qc_col = f"phenotype_qc_{ancestry}"
    lam_col = f"lambda_gc_{ancestry}"
    nca_col = f"n_cases_{ancestry}"

    h2z = pd.to_numeric(m[h2z_col], errors="coerce")
    keep = m["trait_type"].isin(trait_types) & (h2z > min_h2_z)
    if require_qc_pass:
        keep &= m[qc_col] == "PASS"
    if lambda_gc_range is not None:
        lam = pd.to_numeric(m[lam_col], errors="coerce")
        keep &= lam.between(*lambda_gc_range)
    if n_cases_range is not None:
        nca = pd.to_numeric(m[nca_col], errors="coerce")
        keep &= nca.between(*n_cases_range)

    out = m[keep].copy()
    out = out.assign(_h2z=pd.to_numeric(out[h2z_col], errors="coerce"))
    out = out.sort_values("_h2z", ascending=False).drop(columns="_h2z")

    if max_per_category is not None:
        cat = out["category"].fillna("").astype(str)
        exempt = cat.isin(("", "NA"))  # heterogeneous "no category" — don't collapse
        within_rank = out.groupby(cat, sort=False).cumcount()  # 0-based, h2z order
        out = out[exempt | (within_rank < max_per_category)]

    return out.head(n).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Sumstats region access                                                       #
# --------------------------------------------------------------------------- #
def _s3_to_https(path: str) -> str:
    """Map an ``s3://pan-ukb-us-east-1/KEY`` path to its https URL (or pass-through)."""
    prefix = f"s3://{S3_BUCKET}/"
    if path.startswith(prefix):
        return f"{S3_BASE}/{path[len(prefix) :]}"
    return path


def _sumstats_urls(pheno: str | Mapping[str, Any] | pd.Series) -> tuple[str, str]:
    """Resolve a phenotype to ``(data_url, tabix_index_url)`` (https).

    Parameters
    ----------
    pheno
        Either the flat-file ``aws_path`` string, or a manifest row (mapping /
        ``Series``) carrying ``aws_path`` (and optionally ``aws_path_tabix``).

    Returns
    -------
    (str, str)
        The https data URL and its tabix ``.tbi`` index URL. When only the data
        path is known, the index path is derived (Pan-UKBB stores indexes under
        the sibling ``sumstats_flat_files_tabix/`` prefix).
    """
    if isinstance(pheno, str):
        data_path, idx_path = pheno, None
    else:
        data_path = pheno["aws_path"]
        idx_path = pheno.get("aws_path_tabix") if hasattr(pheno, "get") else None
        if idx_path is None and "aws_path_tabix" in pheno:
            idx_path = pheno["aws_path_tabix"]
    data_url = _s3_to_https(data_path)
    if idx_path:
        idx_url = _s3_to_https(idx_path)
    else:
        idx_url = data_url.replace("/sumstats_flat_files/", "/sumstats_flat_files_tabix/") + ".tbi"
    return data_url, idx_url


def _flat_header_columns(data_url: str) -> list[str]:
    """Read a flat file's column header via a small HTTP range request (cached).

    Only the first bgzf block is fetched and decompressed — the header line is
    the first line of the file, so a single ~64 KB range suffices regardless of
    the file's multi-GB total size.
    """
    name = data_url.rsplit("/", 1)[-1]
    cache = CACHE_DIR / "headers" / f"{name}.header.txt"
    if cache.exists():
        return cache.read_text().rstrip("\n").split("\t")
    resp = _session().get(data_url, headers={"Range": "bytes=0-65535"}, timeout=60)
    resp.raise_for_status()
    # Decompress just the first gzip member; the header is well within it.
    with gzip.GzipFile(fileobj=io.BytesIO(resp.content)) as fh:
        first = fh.readline().decode("utf-8").rstrip("\n")
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(first + "\n")
    return first.split("\t")


def _ensure_local_index(idx_url: str) -> Path:
    """Download a remote tabix ``.tbi`` index into the local cache (once).

    Pointing htslib at a *local* index avoids two problems with a remote index:
    it is otherwise re-fetched on every region query (slow for many-window
    workloads), and htslib writes the downloaded index into the current working
    directory (litter). The data file itself is still streamed remotely.
    """
    name = idx_url.rsplit("/", 1)[-1]
    local = CACHE_DIR / "indexes" / name
    if local.exists():
        return local
    resp = _session().get(idx_url, timeout=120)
    resp.raise_for_status()
    local.parent.mkdir(parents=True, exist_ok=True)
    local.write_bytes(resp.content)
    return local


def _tabix_fetch(
    data_url: str, idx_url: str, chrom: str, start: int, end: int, backend: str
) -> list[str]:
    """Return the raw TSV lines for ``chrom:start-end`` via htslib tabix.

    Uses htslib's ``file##idx##index`` addressing so the data file (streamed
    remotely over HTTP range requests) and its tabix index (cached locally, see
    :func:`_ensure_local_index`) can live in different places.
    """
    chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom  # Pan-UKBB is unprefixed
    region = f"{chrom}:{start}-{end}"
    local_index = _ensure_local_index(idx_url)

    if backend == "auto":
        backend = "pysam" if _has_pysam() else "cli"

    if backend == "pysam":
        import pysam  # noqa: PLC0415

        tbx = pysam.TabixFile(data_url, index=str(local_index))
        try:
            return list(tbx.fetch(chrom, start - 1, end))  # pysam is 0-based half-open
        finally:
            tbx.close()

    if backend == "cli":
        if shutil.which("tabix") is None:
            raise RuntimeError(
                "The htslib 'tabix' binary is required for sumstats_region "
                "(default backend). Install htslib, or pass backend='pysam'."
            )
        combined = f"{data_url}##idx##{local_index}"
        proc = subprocess.run(
            ["tabix", combined, region],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"tabix failed for {region}: {proc.stderr.strip()}")
        return [ln for ln in proc.stdout.splitlines() if ln]

    raise ValueError(f"unknown backend: {backend!r} (expected 'auto'|'cli'|'pysam')")


def _has_pysam() -> bool:
    import importlib.util  # noqa: PLC0415

    return importlib.util.find_spec("pysam") is not None


def _parse_sumstats_lines(
    header: list[str],
    lines: list[str],
    ancestries: tuple[str, ...] = ("EUR",),
    *,
    add_z: bool = False,
    add_pval: bool = False,
) -> pd.DataFrame:
    """Parse raw flat-file TSV lines into a tidy per-variant ``DataFrame``.

    Parameters
    ----------
    header
        The flat file's column names (from :func:`_flat_header_columns`).
    lines
        Raw tab-delimited data lines for a region.
    ancestries
        Which ancestries' ``(af, beta, se, neglog10_pval, low_confidence)``
        columns to include. ``"meta"`` selects the meta-analysis track.
    add_z
        Also emit ``z_{anc} = beta_{anc} / se_{anc}``.
    add_pval
        Also emit ``pval_{anc} = 10 ** (-neglog10_pval_{anc})``.

    Returns
    -------
    pandas.DataFrame
        Columns: ``chrom, pos, ref, alt`` then the requested per-ancestry stats.

    Raises
    ------
    KeyError
        If a requested ancestry's effect-size column is absent from ``header``
        (e.g. an ancestry that did not pass QC for this phenotype).
    """
    idx = {c: i for i, c in enumerate(header)}
    for anc in ancestries:
        if f"beta_{anc}" not in idx:
            raise KeyError(f"ancestry {anc!r} not in this flat file (no beta_{anc} column)")

    records = [ln.split("\t") for ln in lines]
    core = {
        "chrom": [r[idx["chr"]] for r in records],
        "pos": [r[idx["pos"]] for r in records],
        "ref": [r[idx["ref"]] for r in records],
        "alt": [r[idx["alt"]] for r in records],
    }
    df = pd.DataFrame(core)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("int64")

    for anc in ancestries:
        for stat in _PER_ANCESTRY_STATS:
            col = f"{stat}_{anc}"
            if col not in idx:
                continue
            values = [r[idx[col]] for r in records]
            if stat == "low_confidence":
                df[col] = pd.Series(values, dtype="object").map(
                    {"true": True, "false": False, "TRUE": True, "FALSE": False}
                )
            else:
                df[col] = pd.to_numeric(pd.Series(values), errors="coerce")
        if add_z:
            df[f"z_{anc}"] = df[f"beta_{anc}"] / df[f"se_{anc}"]
        if add_pval:
            df[f"pval_{anc}"] = 10 ** (-df[f"neglog10_pval_{anc}"])
    return df


def sumstats_region(
    pheno: str | Mapping[str, Any] | pd.Series,
    chrom: str,
    start: int,
    end: int,
    *,
    ancestries: tuple[str, ...] = ("EUR",),
    add_z: bool = True,
    add_pval: bool = False,
    backend: str = "auto",
) -> pd.DataFrame:
    """Fetch dense per-variant summary statistics for a genomic region.

    Streams the region over HTTP range requests via htslib tabix — the phenotype
    flat files are multi-GB but only the requested window is transferred.

    Parameters
    ----------
    pheno
        A flat-file ``aws_path`` string, or a manifest row carrying ``aws_path``.
    chrom, start, end
        Region in **GRCh37/hg19**, 1-based inclusive. ``chrom`` may carry a
        ``chr`` prefix (stripped) — Pan-UKBB contigs are ``1..22, X``.
    ancestries
        Ancestries to return (default EUR). ``"meta"`` = meta-analysis track.
    add_z
        Emit ``z = beta/se`` per ancestry (default True).
    add_pval
        Emit linear ``pval = 10 ** (-neglog10_pval)`` per ancestry.
    backend
        ``"auto"`` (pysam if installed else the ``tabix`` CLI), ``"cli"``, or
        ``"pysam"``.

    Returns
    -------
    pandas.DataFrame
        Tidy per-variant frame (see :func:`_parse_sumstats_lines`).
    """
    data_url, idx_url = _sumstats_urls(pheno)
    header = _flat_header_columns(data_url)
    lines = _tabix_fetch(data_url, idx_url, chrom, start, end, backend)
    return _parse_sumstats_lines(header, lines, ancestries, add_z=add_z, add_pval=add_pval)


def _regions_to_bed(regions: Iterable[tuple[str, int, int]]) -> str:
    """Serialize 1-based-inclusive ``(chrom, start, end)`` regions to BED text.

    BED is 0-based half-open, so ``start`` shifts down by one and ``end`` is kept.
    A leading ``chr`` is stripped (Pan-UKBB contigs are unprefixed).
    """
    rows = []
    for chrom, start, end in regions:
        c = str(chrom)
        c = c[3:] if c.lower().startswith("chr") else c
        rows.append(f"{c}\t{int(start) - 1}\t{int(end)}")
    return ("\n".join(rows) + "\n") if rows else ""


def _tabix_fetch_regions(
    data_url: str, idx_url: str, regions: list[tuple[str, int, int]], backend: str
) -> list[str]:
    """Fetch many regions from one file in a single htslib-tabix pass (``-R bed``)."""
    local_index = _ensure_local_index(idx_url)
    if backend == "auto":
        backend = "pysam" if _has_pysam() else "cli"

    if backend == "pysam":
        import pysam  # noqa: PLC0415

        tbx = pysam.TabixFile(data_url, index=str(local_index))
        try:
            out: list[str] = []
            for chrom, start, end in regions:
                c = str(chrom)
                c = c[3:] if c.lower().startswith("chr") else c
                out.extend(tbx.fetch(c, int(start) - 1, int(end)))  # 0-based half-open
            return out
        finally:
            tbx.close()

    if backend == "cli":
        if shutil.which("tabix") is None:
            raise RuntimeError(
                "The htslib 'tabix' binary is required for sumstats_regions "
                "(default backend). Install htslib, or pass backend='pysam'."
            )
        combined = f"{data_url}##idx##{local_index}"
        fd, bed_path = tempfile.mkstemp(suffix=".bed")
        try:
            with os.fdopen(fd, "w") as fh:
                fh.write(_regions_to_bed(regions))
            proc = subprocess.run(
                ["tabix", "-R", bed_path, combined], capture_output=True, text=True, check=False
            )
        finally:
            os.unlink(bed_path)
        if proc.returncode != 0:
            raise RuntimeError(f"tabix -R failed: {proc.stderr.strip()}")
        return [ln for ln in proc.stdout.splitlines() if ln]

    raise ValueError(f"unknown backend: {backend!r} (expected 'auto'|'cli'|'pysam')")


def sumstats_regions(
    pheno: str | Mapping[str, Any] | pd.Series,
    regions: Iterable[tuple[str, int, int]],
    *,
    ancestries: tuple[str, ...] = ("EUR",),
    add_z: bool = True,
    add_pval: bool = False,
    backend: str = "auto",
) -> pd.DataFrame:
    """Fetch many regions from one phenotype in a single tabix pass.

    The efficient primitive for materializing a training set: one streamed
    ``tabix -R`` call returns every variant across all requested windows, instead
    of one remote query per window. (htslib may merge overlapping regions and
    does not preserve input order; callers assign variants to windows by
    position, so neither matters.)

    Parameters
    ----------
    pheno
        A flat-file ``aws_path`` string, or a manifest row carrying ``aws_path``.
    regions
        Iterable of ``(chrom, start, end)`` in **GRCh37/hg19**, 1-based inclusive.
    ancestries, add_z, add_pval, backend
        As in :func:`sumstats_region`.

    Returns
    -------
    pandas.DataFrame
        Tidy per-variant frame across all regions (see
        :func:`_parse_sumstats_lines`).
    """
    data_url, idx_url = _sumstats_urls(pheno)
    header = _flat_header_columns(data_url)
    lines = _tabix_fetch_regions(data_url, idx_url, list(regions), backend)
    return _parse_sumstats_lines(header, lines, ancestries, add_z=add_z, add_pval=add_pval)
