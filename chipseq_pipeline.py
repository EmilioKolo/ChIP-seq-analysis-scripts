"""
chipseq_pipeline.py
-------------------
Core library for ChIP-seq peak analysis and transcription factor binding site detection.

Pipelines
---------
pipeline_generador()
    Generates binding-site and nearby-gene tables from ChIP-seq BED files and RNA-seq CSV files.

pipeline_meme_chip()
    Generates FASTA sequences from pipeline_generador() outputs for submission to MEME-ChIP.

pipeline_otros_tf()
    Finds binding sites of additional TFs near the primary-TF sites from pipeline_generador().
"""

import logging
import os
from pathlib import Path
from random import shuffle

import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO, motifs
from pyensembl import EnsemblRelease
import urllib.request

try:
    import pyfaidx
except ImportError:
    raise ImportError(
        "pyfaidx is required for automatic genome access. "
        "Install it with: pip install pyfaidx"
    )

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Entrez configuration — set via environment variables or caller code
# ---------------------------------------------------------------------------
Entrez.email = os.environ.get("ENTREZ_EMAIL", "")
Entrez.api_key = os.environ.get("ENTREZ_API_KEY", "")

# ---------------------------------------------------------------------------
# Genome reference constants
# ---------------------------------------------------------------------------
GENOME_DEFAULTS = {
    "human": {"ensembl_release": 75, "genome_name": "hg19", "species": "human"},
    "mouse": {"ensembl_release": 54, "genome_name": "mm9",  "species": "mouse"},
}

# ---------------------------------------------------------------------------
# Genome sequence registry — maps short names to download URLs
# ---------------------------------------------------------------------------

GENOME_CACHE_DIR = Path.home() / ".cache" / "chipseq_pipeline"

# hg19 and hg37 are aliases for GRCh37; both point to the same UCSC file.
GENOME_URLS = {
    "hg19": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
    "hg37": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
    "hg38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    "mm9":  "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz",
    "mm10": "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
}


def get_genome_fasta(genome_name: str, cache_dir: str | Path | None = None) -> "pyfaidx.Fasta":
    """
    Return a pyfaidx.Fasta object for genome_name, downloading and indexing
    the whole-genome FASTA on first use.

    The uncompressed FASTA and its .fai index are stored in cache_dir
    (default: ~/.cache/chipseq_pipeline/). Subsequent calls reuse the cache.

    Parameters
    ----------
    genome_name:
        One of 'hg19', 'hg37' (alias for hg19), 'hg38', 'mm9', 'mm10'.
    cache_dir:
        Override the default cache directory.

    Returns
    -------
    pyfaidx.Fasta
        Random-access FASTA handle. Use genome_fa[chr_n][start:end].seq
        to extract sequences (0-based, half-open).
    """

    key = genome_name.lower()
    if key not in GENOME_URLS:
        raise ValueError(
            f"Genome '{genome_name}' not in registry. "
            f"Available: {list(GENOME_URLS)}"
        )

    cache = Path(cache_dir) if cache_dir else GENOME_CACHE_DIR
    cache.mkdir(parents=True, exist_ok=True)

    # Canonical name — hg37 stores under hg19 since they share the same file
    canonical = "hg19" if key == "hg37" else key
    fasta_path = cache / f"{canonical}.fa"
    gz_path    = cache / f"{canonical}.fa.gz"

    if not fasta_path.exists():
        url = GENOME_URLS[key]
        logger.info("Genome FASTA not found in cache. Downloading %s → %s", url, gz_path)
        logger.info("This is a one-time download (~1 GB compressed). Please be patient.")
        urllib.request.urlretrieve(url, gz_path, reporthook=_download_progress)
        logger.info("Decompressing %s ...", gz_path)
        import gzip, shutil
        with gzip.open(gz_path, "rb") as f_in, open(fasta_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        gz_path.unlink()   # remove the .gz after decompressing
        logger.info("Genome saved to %s", fasta_path)

    # Build .fai index if missing (pyfaidx does this automatically on first open)
    logger.debug("Opening FASTA index: %s", fasta_path)
    return pyfaidx.Fasta(str(fasta_path), as_raw=True)


def _download_progress(block_num: int, block_size: int, total_size: int) -> None:
    """urllib reporthook that logs download progress every ~5%."""
    if total_size <= 0:
        return
    pct = block_num * block_size * 100 / total_size
    if block_num == 0 or int(pct) % 5 == 0:
        logger.info("Download progress: %.1f%%", min(pct, 100))


def load_genome(organism: str, ensembl_release: int | None = None) -> tuple[EnsemblRelease, str]:
    """
    Return an EnsemblRelease object and genome name string for the given organism.

    Parameters
    ----------
    organism:
        Either 'human' or 'mouse'.
    ensembl_release:
        Ensembl release number. If None, the default for the organism is used
        (75 for human / hg19, 54 for mouse / mm9).

    Returns
    -------
    genome : EnsemblRelease
    genome_name : str
        Short genome name (e.g. 'hg19', 'mm9'), or '<species>_r<release>'
        when overriding the default release.
    """
    key = organism.lower()
    if key not in GENOME_DEFAULTS:
        raise ValueError(f"Organism '{organism}' not recognised. Choose 'human' or 'mouse'.")

    defaults = GENOME_DEFAULTS[key]
    release = ensembl_release if ensembl_release is not None else defaults["ensembl_release"]
    genome_name = defaults["genome_name"] if release == defaults["ensembl_release"] else f"{defaults['species']}_r{release}"

    logger.info("Loading EnsemblRelease %d for %s (genome name: %s).", release, organism, genome_name)
    genome = EnsemblRelease(release, species=defaults["species"])
    return genome, genome_name


# ===========================================================================
# Top-level pipeline functions
# ===========================================================================


def pipeline_generador(
    bed_file: str,
    rnaseq_file: str,
    output_prefix: str,
    genome: EnsemblRelease,
    genome_name: str,
    pssm,
    score_cutoff: float,
    genome_cache_dir: str | Path | None = None,
    path_bed: str = "",
    path_out: str = "",
    path_fasta: str = "",
    dist_max: int = 100_000,
    binding_sites: list[str] | None = None,
    test_mode: int = 0,
    id_col_rnaseq: int = 0,
    updown_col_rnaseq: int = 3,
    rnaseq_sep: str = ";",
    translate_col: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Main generator pipeline: BED peaks + RNA-seq → binding-site and gene tables.

    Parameters
    ----------
    bed_file:
        Name of the BED file (without extension) located in path_bed.
    rnaseq_file:
        Name of the RNA-seq CSV file (without extension) located in path_bed.
    output_prefix:
        Base name for all output files.
    genome:
        EnsemblRelease object used to look up nearby genes.
    genome_name:
        Short genome identifier used to load per-chromosome FASTA files.
    pssm:
        Biopython PSSM object for the primary transcription factor.
    score_cutoff:
        Minimum PSSM score to call a binding site.
    path_bed:
        Directory containing the BED and RNA-seq input files.
    path_out:
        Directory where output CSVs are written.
    path_fasta:
        Directory containing per-chromosome FASTA files.
    dist_max:
        Maximum distance (bp) from a peak centre to search for nearby genes.
    binding_sites:
        List of exact sequences treated as confirmed binding sites (list-based search).
        Defaults to an empty list (no list-based search).
    test_mode:
        If > 0, run on a random subset of this many peaks (for debugging).
    id_col_rnaseq:
        Column index for the gene identifier in the RNA-seq file.
    updown_col_rnaseq:
        Column index for the log(fold-change) value in the RNA-seq file.
    rnaseq_sep:
        Field separator used in the RNA-seq CSV file.
    translate_col:
        If not None, the column index containing gene *names* to translate to
        Ensembl IDs using *genome*. When None, identifiers are used as-is.

    Returns
    -------
    df_peaks, df_binding_sites, df_genes : pd.DataFrame
    """
    if binding_sites is None:
        binding_sites = []

    translate = [translate_col, genome] if translate_col is not None else []

    # ------------------------------------------------------------------
    # Load inputs
    # ------------------------------------------------------------------
    df_peaks = load_bed(bed_file, path_dir=path_bed)
    dict_rnaseq = load_rnaseq(
        rnaseq_file,
        path_dir=path_bed,
        sep=rnaseq_sep,
        id_col=id_col_rnaseq,
        updown_col=updown_col_rnaseq,
        translate=translate,
    )

    if test_mode > 0:
        logger.info("Test mode: shuffling and selecting %d peaks.", test_mode)
        rows = df_peaks.to_dict("records")
        shuffle(rows)
        df_peaks = pd.DataFrame(rows[:test_mode])

    peaks = df_peaks.to_dict("records")

    # ------------------------------------------------------------------
    # Load FASTA sequences for all chromosomes present in the BED file
    # ------------------------------------------------------------------
    dict_chr_seq = load_chr_sequences(
        peaks, genome_name, path_fasta=path_fasta, genome_cache_dir=genome_cache_dir
    )

    # ------------------------------------------------------------------
    # Process each peak
    # ------------------------------------------------------------------
    rows_peaks = []
    rows_su = []
    rows_genes = []

    for idx, peak in enumerate(peaks):
        chr_n = peak["chr_n"]
        pos_ini = int(peak["pos_ini"])
        pos_end = int(peak["pos_end"])

        if chr_n not in dict_chr_seq:
            logger.warning("Chromosome %s not found in FASTA dict; skipping peak %d.", chr_n, idx)
            continue

        # Nearby genes
        genes_near = find_nearby_genes(chr_n, pos_ini, pos_end, dist_max, genome)
        genes_near, n_total, n_up, n_down = add_rnaseq_info(genes_near, dict_rnaseq)

        # Sequence and binding sites
        seq_peak = extract_sequence(dict_chr_seq[chr_n], pos_ini, pos_end)
        su_list = find_binding_sites_list(seq_peak, chr_n, binding_sites, pos_ini - 1)
        su_pssm = find_binding_sites_pssm(seq_peak, chr_n, pssm, score_cutoff, pos_ini - 1)

        # Annotate binding sites with gene context
        for su in su_list:
            su.update({
                "score_pssm": None,
                "n_genes_near": len(genes_near),
                "n_updown_total": int(n_total),
                "n_upreg": int(n_up),
                "n_downreg": int(n_down),
            })
        for su in su_pssm:
            su.update({
                "n_genes_near": len(genes_near),
                "n_updown_total": int(n_total),
                "n_upreg": int(n_up),
                "n_downreg": int(n_down),
            })

        # Annotate peak row
        peak.update({
            "n_genes": len(genes_near),
            "n_su": len(su_list) + len(su_pssm),
            "n_su_list": len(su_list),
            "n_su_pssm": len(su_pssm),
            "n_updown_total": int(n_total),
            "n_upreg": int(n_up),
            "n_downreg": int(n_down),
        })

        rows_peaks.append(peak)
        rows_su.extend(su_list + su_pssm)
        rows_genes.extend(genes_near)

        if (idx + 1) % 100 == 0 or idx == 0:
            logger.info("Progress: %d / %d peaks processed.", idx + 1, len(peaks))

    # ------------------------------------------------------------------
    # Deduplicate genes and build output DataFrames
    # ------------------------------------------------------------------
    dict_chr_seq.clear()

    df_genes_all = pd.DataFrame(rows_genes)
    if not df_genes_all.empty:
        df_genes_all = df_genes_all.drop_duplicates(subset=["gene_id"])

    df_su = pd.DataFrame(rows_su) if rows_su else pd.DataFrame(
        columns=["chr_n", "pos_ini", "pos_end", "seq", "source", "score_pssm",
                 "n_genes_near", "n_updown_total", "n_upreg", "n_downreg"]
    )
    df_peaks_out = pd.DataFrame(rows_peaks) if rows_peaks else pd.DataFrame(
        columns=["chr_n", "pos_ini", "pos_end", "n_genes", "n_su", "n_su_list",
                 "n_su_pssm", "n_updown_total", "n_upreg", "n_downreg"]
    )

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    save_dataframe(df_su, f"{output_prefix}_binding_sites", path_out=path_out)
    save_dataframe(df_genes_all, f"{output_prefix}_genes", path_out=path_out)
    save_dataframe(df_peaks_out, f"{output_prefix}_peaks", path_out=path_out)

    return df_peaks_out, df_su, df_genes_all


def pipeline_meme_chip(
    sites_file: str,
    output_name: str,
    genome_name: str,
    genome_cache_dir: str | Path | None = None,
    site_length: int = 0,
    col_sites: list[int] | None = None,
    filter_rules: list[tuple] | None = None,
    default_pass_filter: bool = True,
    path_sites: str = "",
    path_out: str = "",
    path_fasta: str = "",
) -> pd.DataFrame:
    """
    Generate a FASTA file of peak/binding-site sequences for MEME-ChIP input.

    Parameters
    ----------
    sites_file:
        Name of the input CSV file (without extension) with peak or binding-site coordinates.
    output_name:
        Name of the output FASTA file (without extension).
    genome_name:
        Short genome name used to locate per-chromosome FASTA files in path_fasta.
    site_length:
        If > 0, all sequences are resized to this length (centred on the original site).
        If 0, sequences are returned at their original length.
    col_sites:
        Column indices for [chr_n, pos_ini, pos_end] in the input CSV.
        Defaults to [0, 1, 2].
    filter_rules:
        List of (column_index, expected_value) tuples used to include/exclude sites.
    default_pass_filter:
        If True, all sites pass by default and filter_rules *exclude* sites.
        If False, no site passes by default and filter_rules *include* sites.
    path_sites:
        Directory containing the input CSV.
    path_out:
        Directory where the output FASTA is written.
    path_fasta:
        Directory containing per-chromosome FASTA files.

    Returns
    -------
    df_filtered : pd.DataFrame
        Table of sites that passed the filter (with coordinates as loaded/resized).
    """
    if col_sites is None:
        col_sites = [0, 1, 2]
    if filter_rules is None:
        filter_rules = []

    df_sites = load_csv(sites_file, path_dir=path_sites)

    # ------------------------------------------------------------------
    # Filter sites
    # ------------------------------------------------------------------
    logger.info("Filtering %d sites.", len(df_sites))
    mask = df_sites.apply(lambda row: _apply_filter(row.tolist(), filter_rules, default_pass_filter), axis=1)
    df_filtered = df_sites[mask].copy()
    df_filtered = df_filtered.iloc[:, col_sites].copy()
    df_filtered.columns = ["chr_n", "pos_ini", "pos_end"]
    df_filtered = df_filtered[df_filtered["chr_n"] != ""]

    # ------------------------------------------------------------------
    # Optionally resize to uniform length
    # ------------------------------------------------------------------
    if site_length > 0:
        logger.info("Resizing sites to %d bp.", site_length)
        df_filtered[["pos_ini", "pos_end"]] = df_filtered.apply(
            lambda row: pd.Series(resize_site(int(row["pos_ini"]), int(row["pos_end"]), site_length)),
            axis=1,
        )

    # ------------------------------------------------------------------
    # Extract sequences
    # ------------------------------------------------------------------
    logger.info("Extracting sequences for %d sites.", len(df_filtered))
    sites = df_filtered.to_dict("records")
    dict_chr_seq = load_chr_sequences(
        sites, genome_name, path_fasta=path_fasta, genome_cache_dir=genome_cache_dir
    )

    fasta_records: dict[str, str] = {}
    for i, site in enumerate(sites):
        chr_n = site["chr_n"]
        if chr_n not in dict_chr_seq:
            logger.warning("Chromosome %s not found; skipping site %d.", chr_n, i)
            continue
        seq = extract_sequence(dict_chr_seq[chr_n], int(site["pos_ini"]), int(site["pos_end"]))
        seq_id = f"{i + 1}_{chr_n}_{site['pos_ini']}_{site['pos_end']}"
        fasta_records[seq_id] = str(seq)
        if (i + 1) % 1000 == 0 or i == 0:
            logger.info("Progress: %d / %d sites.", i + 1, len(sites))

    dict_chr_seq.clear()
    save_fasta(fasta_records, output_name, path_out=path_out)

    return df_filtered


def pipeline_otros_tf(
    sites_file: str,
    genes_file: str,
    output_prefix: str,
    genome_name: str,
    pwm_files: list[str],
    genome_cache_dir: str | Path | None = None,
    dist_sites: int = 250,
    dist_max_genes: int = 100_000,
    pwm_names: list[str] | None = None,
    path_sites: str = "",
    path_genes: str = "",
    path_out: str = "",
    path_fasta: str = "",
    path_pwm: str = "",
    score_mult: float = 0.9,
    score_report: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Find binding sites of additional TFs near the primary-TF sites.

    Parameters
    ----------
    sites_file:
        Name of the primary-TF sites CSV (without extension) from pipeline_generador().
    genes_file:
        Name of the nearby-genes CSV (without extension) from pipeline_generador().
    output_prefix:
        Base name for all output files.
    genome_name:
        Short genome name for FASTA loading.
    pwm_files:
        List of PWM/PCM filenames (without path) for the additional TFs.
    dist_sites:
        Distance (bp) to expand each primary-TF site when searching for other TFs.
        Use 0 or negative to search only within the site itself.
    dist_max_genes:
        Maximum distance to the gene +1 when reporting nearby genes.
    pwm_names:
        Display names for the additional TFs (same order as pwm_files).
    path_sites / path_genes / path_out / path_fasta / path_pwm:
        Directories for input files, output files, FASTA files, and PWM files.
    score_mult:
        Fraction of PSSM maximum score used as the binding-site score cutoff.
    score_report:
        If True, print max-score information for every PSSM (useful for cutoff calibration).

    Returns
    -------
    df_sites_out : pd.DataFrame
        Primary-TF sites annotated with nearby genes and co-occupancy flags.
    df_other_tf : pd.DataFrame
        All binding sites found for the additional TFs, with primary-TF site context.
    """
    if pwm_names is None:
        pwm_names = []

    # ------------------------------------------------------------------
    # Load PSSMs
    # ------------------------------------------------------------------
    pssm_list = []
    for pwm_file in pwm_files:
        pssm = load_pssm(pwm_file, path_dir=path_pwm)
        pssm_list.append(pssm)

    if score_report:
        for name, pssm in zip(pwm_names or pwm_files, pssm_list):
            logger.info("PSSM score report — %s: max=%.4f, cutoff@%.0f%%=%.4f",
                        name, pssm.max, score_mult * 100, pssm.max * score_mult)

    # ------------------------------------------------------------------
    # Load inputs
    # ------------------------------------------------------------------
    df_sites, header_sites = load_csv(sites_file, path_dir=path_sites, return_header=True)
    dict_genes = load_dict_csv(genes_file, key_col=1, path_dir=path_genes)

    logger.info("Processing %d primary-TF sites.", len(df_sites))
    sites = df_sites.to_dict("records")
    raw_cols = list(df_sites.columns)

    dict_chr_seq = load_chr_sequences(
        sites, genome_name, path_fasta=path_fasta, genome_cache_dir=genome_cache_dir
    )

    rows_sites_out = []
    rows_other_tf = []
    skipped = 0

    for i, site in enumerate(sites):
        chr_n = site.get("chr_n", site[raw_cols[0]])
        pos_ini_orig = int(site.get("pos_ini", site[raw_cols[1]]))
        pos_end_orig = int(site.get("pos_end", site[raw_cols[2]]))

        if chr_n not in dict_chr_seq:
            logger.warning("Chromosome %s not in FASTA dict; skipping site %d.", chr_n, i)
            skipped += 1
            continue

        pos_ini_search = pos_ini_orig - dist_sites if dist_sites > 0 else pos_ini_orig
        pos_end_search = pos_end_orig + dist_sites if dist_sites > 0 else pos_end_orig

        seq_range = extract_sequence(dict_chr_seq[chr_n], pos_ini_search, pos_end_search)

        other_tf_hits, genes_near = _find_other_tf_sites(
            chr_n, pos_ini_orig, pos_end_orig,
            dict_genes, dist_max_genes,
            seq_range, pssm_list,
            pssm_names=pwm_names,
            score_mult=score_mult,
            pos_ini_search=pos_ini_search,
        )

        genes_str = ",".join(f"{g['gene_id']}_{g['updown']}" for g in genes_near)

        site_row = dict(site)
        site_row["str_genes_rnaseq_near"] = genes_str
        tf_flags = {name: 0 for name in (pwm_names or pwm_files)}
        for hit in other_tf_hits:
            tf_name = hit.get("tf_name", "")
            if tf_name in tf_flags:
                tf_flags[tf_name] = 1
        site_row.update(tf_flags)
        rows_sites_out.append(site_row)

        for hit in other_tf_hits:
            combined = dict(hit)
            combined.update(site_row)
            rows_other_tf.append(combined)

        if (i + 1) % 250 == 0 or i == 0:
            logger.info("Progress: %d / %d sites.", i + 1, len(sites))

    dict_chr_seq.clear()
    logger.info("Skipped %d sites (chromosome not found in FASTA).", skipped)

    df_sites_out = pd.DataFrame(rows_sites_out)
    df_other_tf = pd.DataFrame(rows_other_tf)

    save_dataframe(df_sites_out, f"{output_prefix}_sites_genes", path_out=path_out)
    save_dataframe(df_other_tf, f"{output_prefix}_other_tf", path_out=path_out)

    return df_sites_out, df_other_tf


# ===========================================================================
# Core biological functions
# ===========================================================================


def find_nearby_genes(
    chr_n: str,
    pos_ini: int,
    pos_end: int,
    dist_max: int,
    genome: EnsemblRelease,
) -> list[dict]:
    """
    Return protein-coding genes within dist_max bp of the peak defined by
    chr_n, pos_ini, pos_end.

    Each returned gene dict has keys: gene_id, chr_n, pos0, biotype.
    """
    contig = _chr_to_contig(chr_n)
    try:
        genes_raw = genome.genes_at_locus(contig, pos_ini - dist_max, pos_end + dist_max)
    except Exception:
        logger.exception("Error querying genes at locus %s:%d-%d.", chr_n, pos_ini, pos_end)
        return []

    result = []
    for gene in genes_raw:
        if gene.biotype != "protein_coding":
            continue
        pos0 = gene.start if gene.strand == "+" else gene.end
        result.append({
            "gene_id": gene.gene_id,
            "chr_n": chr_n,
            "pos0": pos0,
            "biotype": gene.biotype,
        })
    return result


def add_rnaseq_info(
    genes: list[dict],
    dict_rnaseq: dict[str, float],
) -> tuple[list[dict], int, int, int]:
    """
    Annotate each gene dict with its RNA-seq fold-change value.

    Returns the annotated list plus counts: (n_total, n_up, n_down).
    n_total counts genes with a non-zero fold-change in the RNA-seq data.
    """
    n_total = n_up = n_down = 0
    result = []
    for gene in genes:
        gene_out = dict(gene)
        gene_id = gene["gene_id"]
        if gene_id in dict_rnaseq:
            fc = dict_rnaseq[gene_id]
            gene_out["fold_change"] = fc
            n_total += 1
            if fc > 0:
                n_up += 1
            elif fc < 0:
                n_down += 1
            else:
                logger.warning("Gene %s has fold_change=0.", gene_id)
        else:
            gene_out["fold_change"] = None
        result.append(gene_out)
    return result, n_total, n_up, n_down


def find_binding_sites_list(
    seq_peak: str,
    chr_n: str,
    site_sequences: list[str],
    pos_ini_ref: int,
) -> list[dict]:
    """
    Search seq_peak for exact matches to each sequence in site_sequences,
    including the reverse complement of each.

    Returns a list of hit dicts: chr_n, pos_ini, pos_end, seq, source='list'.
    """
    hits = []
    for seq in site_sequences:
        hits.extend(_find_one_binding_site(seq_peak, seq, chr_n, pos_ini_ref))
    return hits


def find_binding_sites_pssm(
    seq_peak,
    chr_n: str,
    pssm,
    score_cutoff: float,
    pos_ini_ref: int,
) -> list[dict]:
    """
    Scan seq_peak with a Biopython PSSM and return all hits above score_cutoff.

    Returns a list of hit dicts: chr_n, pos_ini, pos_end, seq, source='pssm', score_pssm.
    """
    hits = []
    pssm_len = len(pssm.consensus)
    if len(seq_peak) < pssm_len:
        return hits

    try:
        for position, score in pssm.search(seq_peak, threshold=score_cutoff):
            seq_hit = seq_peak[position : position + pssm_len]
            if position < 0:
                pos_ini = pos_ini_ref + len(seq_peak) + position
            else:
                pos_ini = position + pos_ini_ref
            hits.append({
                "chr_n": chr_n,
                "pos_ini": pos_ini,
                "pos_end": pos_ini + pssm_len - 1,
                "seq": str(seq_hit),
                "source": "pssm",
                "score_pssm": float(score),
            })
    except Exception:
        logger.exception("Error scanning PSSM on sequence starting at ref %d.", pos_ini_ref)

    return hits


def extract_sequence(chr_seq, pos_ini: int, pos_end: int) -> str:
    """
    Return the subsequence for [pos_ini, pos_end] (1-based, inclusive).

    Accepts both Biopython Seq objects (legacy path) and pyfaidx sequence
    objects (auto-genome path) — the slice arithmetic is the same since
    pyfaidx.as_raw=True returns plain strings via 0-based slicing.
    """
    return chr_seq[pos_ini - 1 : pos_end]


def resize_site(pos_ini: int, pos_end: int, target_length: int) -> tuple[int, int]:
    """
    Expand or contract [pos_ini, pos_end] symmetrically to target_length bp.
    Ties are broken by extending pos_end first.
    """
    current = pos_end - pos_ini
    lo, hi = pos_ini, pos_end

    if target_length > current:
        parity = False
        while target_length > (hi - lo):
            if parity:
                hi += 1
            else:
                lo -= 1
            parity = not parity
    elif target_length < current:
        parity = False
        while target_length < (hi - lo):
            if parity:
                lo += 1
            else:
                hi -= 1
            parity = not parity

    if target_length != (hi - lo):
        logger.error("resize_site produced unexpected result: target=%d, got=%d.", target_length, hi - lo)

    return lo, hi


# ===========================================================================
# I/O functions
# ===========================================================================


def load_bed(
    name: str,
    path_dir: str = "",
    sep: str = "\t",
    ext: str = ".bed",
    n_cols: int = 3,
) -> pd.DataFrame:
    """
    Load a BED file and return a DataFrame with columns chr_n, pos_ini, pos_end
    (plus any additional columns up to n_cols).
    """
    filepath = _build_path(path_dir, name, ext)
    logger.info("Loading BED file: %s", filepath)
    col_names = ["chr_n", "pos_ini", "pos_end"] + [f"col_{i}" for i in range(3, n_cols)]
    try:
        df = pd.read_csv(filepath, sep=sep, header=None, usecols=range(n_cols), names=col_names)
    except FileNotFoundError:
        raise FileNotFoundError(f"BED file not found: {filepath}")
    return df


def load_csv(
    name: str,
    path_dir: str = "",
    sep: str = ";",
    ext: str = ".csv",
    return_header: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, list[str]]:
    """
    Load a semicolon-separated CSV and return a DataFrame.
    If return_header is True, also return the list of column names.
    """
    filepath = _build_path(path_dir, name, ext)
    logger.info("Loading CSV: %s", filepath)
    try:
        df = pd.read_csv(filepath, sep=sep)
    except FileNotFoundError:
        raise FileNotFoundError(f"CSV file not found: {filepath}")
    if return_header:
        return df, list(df.columns)
    return df


def load_dict_csv(
    name: str,
    key_col: int = 0,
    path_dir: str = "",
    sep: str = ";",
    ext: str = ".csv",
) -> dict[str, list[list]]:
    """
    Load a CSV and return a dict keyed by the values in key_col.
    Each value is a list of rows (as lists) sharing that key.
    """
    filepath = _build_path(path_dir, name, ext)
    logger.info("Loading dict CSV: %s (key_col=%d)", filepath, key_col)
    try:
        df = pd.read_csv(filepath, sep=sep, header=0)
    except FileNotFoundError:
        raise FileNotFoundError(f"CSV file not found: {filepath}")

    result: dict[str, list] = {}
    for row in df.itertuples(index=False):
        row_list = list(row)
        key = str(row_list[key_col])
        result.setdefault(key, []).append(row_list)
    return result


def load_rnaseq(
    name: str,
    path_dir: str = "",
    sep: str = ";",
    ext: str = ".csv",
    id_col: int = 0,
    updown_col: int = 3,
    translate: list | None = None,
) -> dict[str, float]:
    """
    Load an RNA-seq results file and return a dict mapping gene ID → fold-change.

    The file is expected to contain log2(fold-change) values; these are converted
    to signed fold-changes with sign(fc) * 2**|fc|.

    Parameters
    ----------
    translate:
        If non-empty, should be [gene_name_col_index, EnsemblRelease_object].
        Gene names in that column will be translated to Ensembl IDs.
    """
    if translate is None:
        translate = []

    filepath = _build_path(path_dir, name, ext)
    logger.info("Loading RNA-seq file: %s", filepath)

    result: dict[str, float] = {}
    found = total = 0

    try:
        df = pd.read_csv(filepath, sep=sep, header=0)
    except FileNotFoundError:
        raise FileNotFoundError(f"RNA-seq file not found: {filepath}")

    for _, row in df.iterrows():
        row_list = list(row)
        raw_id = str(row_list[id_col])
        try:
            fc_abs = _fold_change_from_log(float(row_list[updown_col]))
        except (ValueError, TypeError):
            logger.error("Could not parse fold-change for row: %s", row_list)
            continue

        for gene_id in raw_id.split("|"):
            gene_id = gene_id.strip()
            if not gene_id:
                continue

            if not translate:
                _insert_fc(result, gene_id, fc_abs)
                found += 1
            else:
                gene_name_col, genome = translate[0], translate[1]
                gene_name = str(row_list[gene_name_col])
                try:
                    ensembl_ids = genome.gene_ids_of_gene_name(gene_name)
                    for eid in ensembl_ids:
                        _insert_fc(result, eid, fc_abs)
                        found += 1
                except Exception:
                    logger.warning("No Ensembl ID found for gene name '%s'.", gene_name)

        total += 1

    logger.info("RNA-seq: %d total rows, %d identifiers mapped, %d unique genes.", total, found, len(result))
    return result


def load_pssm(name: str, path_dir: str = "", pseudocounts: float = 0.5, solo_pssm: bool = True):
    """
    Parse a HOCOMOCO PCM file and return a Biopython PSSM object.

    Parameters
    ----------
    solo_pssm:
        If True (default), return only the PSSM.
        If False, return (motif, pwm, pssm).
    """
    filepath = _build_path(path_dir, name, "")
    logger.info("Loading PSSM: %s", filepath)
    try:
        with open(filepath) as fh:
            record = motifs.parse(fh, "pfm-four-columns")
    except FileNotFoundError:
        raise FileNotFoundError(f"PWM file not found: {filepath}")

    if len(record) == 0:
        raise ValueError(f"No motifs found in {filepath}.")
    if len(record) > 1:
        logger.warning("%s contains %d motifs; using only the first.", name, len(record))

    m = record[0]
    pwm = m.counts.normalize(pseudocounts=pseudocounts)
    pssm = pwm.log_odds()

    return pssm if solo_pssm else (m, pwm, pssm)


def load_chr_sequences(
    sites: list[dict],
    genome_name: str,
    path_fasta: str = "",
    chr_key: str = "chr_n",
    genome_cache_dir: str | Path | None = None,
) -> dict[str, object]:
    """
    Return a dict mapping chr_n → sequence accessor for every chromosome in sites.

    If path_fasta is provided, the old behaviour is preserved: one .fasta file
    per chromosome is expected at path_fasta/<genome_name>_<chr_n>.fasta.

    If path_fasta is empty, the genome is loaded via get_genome_fasta(), which
    downloads and indexes the whole-genome FASTA automatically on first use.
    Supported genome names: hg19, hg37 (=hg19), hg38, mm9, mm10.
    """
    if path_fasta:
        # Legacy mode — per-chromosome FASTA files (original behaviour)
        chr_set = {str(s[chr_key]) for s in sites}
        result = {}
        for chr_n in sorted(chr_set):
            fasta_path = Path(path_fasta) / f"{genome_name}_{chr_n}.fasta"
            try:
                record = SeqIO.read(str(fasta_path), "fasta")
                result[chr_n] = record.seq
                logger.debug("Loaded FASTA: %s", fasta_path)
            except FileNotFoundError:
                logger.error("FASTA file not found: %s", fasta_path)
            except Exception:
                logger.exception("Error reading FASTA file: %s", fasta_path)
        return result

    # Auto mode — whole-genome FASTA via pyfaidx
    genome_fa = get_genome_fasta(genome_name, cache_dir=genome_cache_dir)

    chr_set = {str(s[chr_key]) for s in sites}
    result = {}
    for chr_n in sorted(chr_set):
        # pyfaidx keys are the FASTA record names; UCSC uses 'chr1', 'chrX', etc.
        if chr_n in genome_fa:
            result[chr_n] = genome_fa[chr_n]
        else:
            logger.warning("Chromosome '%s' not found in genome FASTA.", chr_n)
    return result


def save_dataframe(df: pd.DataFrame, name: str, path_out: str = "", sep: str = ";", ext: str = ".csv") -> None:
    """Save a DataFrame as a semicolon-separated CSV."""
    filepath = _build_path(path_out, name, ext)
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(filepath, sep=sep, index=False)
    logger.info("Saved: %s (%d rows)", filepath, len(df))


def save_fasta(records: dict[str, str], name: str, path_out: str = "", ext: str = ".fasta") -> None:
    """Save a dict of {sequence_id: sequence} as a FASTA file."""
    filepath = _build_path(path_out, name, ext)
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as fh:
        for seq_id, seq in records.items():
            fh.write(f">{seq_id}\n{seq}\n")
    logger.info("Saved FASTA: %s (%d sequences)", filepath, len(records))


# ===========================================================================
# Private helpers
# ===========================================================================


def _build_path(directory: str, name: str, ext: str) -> str:
    """Construct a file path from directory, base name, and extension."""
    filename = name + ext if ext and not name.endswith(ext) else name
    if directory:
        return str(Path(directory) / filename)
    return filename


def _chr_to_contig(chr_n: str) -> str:
    """Convert 'chr1' → '1', 'chrX' → 'X', etc."""
    return chr_n[3:] if chr_n.startswith("chr") else chr_n


def _fold_change_from_log(log_fc: float, base: float = 2.0) -> float:
    """Convert log2(fold-change) to signed fold-change: sign(x) * base^|x|."""
    return float(np.sign(log_fc)) * (base ** abs(log_fc))


def _insert_fc(d: dict, gene_id: str, fc: float) -> None:
    """Insert fold-change into dict, warning on conflicting duplicate entries."""
    if gene_id not in d:
        d[gene_id] = fc
    elif d[gene_id] != fc:
        logger.error(
            "Duplicate gene ID '%s' with conflicting fold-changes: existing=%.4f, new=%.4f.",
            gene_id, d[gene_id], fc,
        )


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence string."""
    from Bio.Seq import Seq
    return str(Seq(seq).reverse_complement())


def _find_one_binding_site(
    seq_reference: str,
    seq_query: str,
    chr_n: str,
    pos_ini_ref: int,
) -> list[dict]:
    """
    Find all occurrences of seq_query or its reverse complement in seq_reference.
    Returns hit dicts: chr_n, pos_ini, pos_end, seq, source='list'.
    """
    hits = []
    query_upper = seq_query.upper()
    revcomp_upper = _reverse_complement(seq_query).upper()
    seq_len = len(seq_query)

    for i in range(len(seq_reference) - seq_len + 1):
        window = str(seq_reference[i : i + seq_len]).upper()
        if window == query_upper or window == revcomp_upper:
            hits.append({
                "chr_n": str(chr_n),
                "pos_ini": pos_ini_ref + i + 1,
                "pos_end": pos_ini_ref + i + seq_len,
                "seq": window,
                "source": "list",
                "score_pssm": None,
            })
    return hits


def _apply_filter(row: list, rules: list[tuple], default_pass: bool) -> bool:
    """
    Return True if row passes the filter defined by rules.

    Each rule is (column_index, expected_value). If default_pass is True,
    a matching rule *excludes* the row; if False, a matching rule *includes* it.
    """
    result = default_pass
    for col_idx, expected in rules:
        if row[col_idx] == expected:
            result = not default_pass
    return result


def _find_other_tf_sites(
    chr_n: str,
    pos_ini: int,
    pos_end: int,
    dict_genes: dict,
    dist_genes: int,
    seq_range,
    pssm_list: list,
    pssm_names: list[str],
    score_mult: float,
    pos_ini_search: int,
) -> tuple[list[dict], list[dict]]:
    """
    Find binding sites for a list of PSSMs and nearby genes for a single primary-TF site.
    Returns (other_tf_hits, genes_near).
    """
    genes_near = []
    if chr_n in dict_genes:
        for gene_row in dict_genes[chr_n]:
            # gene_row format: [gene_id, chr_n, pos0, biotype, fold_change, ...]
            if len(gene_row) > 4 and gene_row[4] != "" and gene_row[4] is not None:
                try:
                    pos0 = int(gene_row[2])
                    fc = float(gene_row[4])
                    updown = "upreg" if fc < 0 else "downreg"
                    if (pos_ini - dist_genes) < pos0 < (pos_end + dist_genes):
                        genes_near.append({"gene_id": gene_row[0], "updown": updown})
                except (ValueError, TypeError):
                    logger.warning("Could not parse gene row: %s", gene_row)

    other_tf_hits = []
    for i, pssm in enumerate(pssm_list):
        name = pssm_names[i] if i < len(pssm_names) else str(i)
        cutoff = score_mult * pssm.max
        hits = find_binding_sites_pssm(seq_range, chr_n, pssm, cutoff, pos_ini_search - 1)
        for hit in hits:
            hit["tf_name"] = name
        other_tf_hits.extend(hits)

    return other_tf_hits, genes_near