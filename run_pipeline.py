#!/usr/bin/env python3
"""
run_pipeline.py
---------------
Command-line entry point for the ChIP-seq binding-site analysis pipeline.

Subcommands
-----------
  generador     Load BED + RNA-seq → binding sites and nearby genes.
  memechip      Generate FASTA sequences for MEME-ChIP input.
  otros-tf      Find co-occupancy of additional TFs near primary-TF sites.

Examples
--------
# Run the full generator pipeline (human, NKX2-5 defaults):
python run_pipeline.py generador \\
    --bed-file Anderson2018 \\
    --rnaseq-file Anderson2018_RNAseq_source \\
    --output-prefix anderson_full \\
    --organism human \\
    --path-bed /data/input \\
    --path-out /data/output/100kpb \\
    --path-fasta /data/genomes \\
    --path-pwm /data/PWM_human \\
    --pssm-file NKX25_HUMAN.H11MO.0.B.pcm \\
    --binding-sites GCAAGTG GGAAGTG GAAAGTG \\
    --dist-max 100000

# Run MEME-ChIP FASTA generation:
python run_pipeline.py memechip \\
    --sites-file anderson_full_peaks \\
    --output-name anderson_full_peaks_fasta \\
    --organism human \\
    --path-sites /data/output/100kpb \\
    --path-out /data/output/100kpb \\
    --path-fasta /data/genomes \\
    --site-length 500

# Run co-occupancy analysis for additional TFs:
python run_pipeline.py otros-tf \\
    --sites-file anderson_full_binding_sites \\
    --genes-file anderson_full_genes \\
    --output-prefix anderson_full_binding_sites_other_tf \\
    --organism human \\
    --pwm-files NKX25_HUMAN.H11MO.0.B.pcm TBX20_HUMAN.H11MO.0.D.pcm \\
    --pwm-names NKX25 TBX20 \\
    --path-sites /data/output/100kpb \\
    --path-genes /data/output/100kpb \\
    --path-out /data/output \\
    --path-fasta /data/genomes \\
    --path-pwm /data/PWM_human \\
    --dist-sites 250
"""

import argparse
import logging
import sys

import chipseq_pipeline as pipeline

# ---------------------------------------------------------------------------
# Default NKX2-5 PSSM filenames (override with --pssm-file)
# ---------------------------------------------------------------------------
DEFAULT_PSSM = {
    "human": "NKX25_HUMAN.H11MO.0.B.pcm",
    "mouse": "NKX25_MOUSE.H11MO.0.A.pcm",
}

# Default confirmed NKX2-5 binding site sequences
DEFAULT_BINDING_SITES = [
    "GCAAGTG", "GGAAGTG", "GAAAGTG", "ATAAGTG",
    "GTAAGTG", "CTAAGTG", "TCAAGTG", "TGAAGTG",
    "TAAAGTG", "TTAAGTG",
]

# Default additional TF PWM files and display names per organism
DEFAULT_PWM_HUMAN = [
    "NKX25_HUMAN.H11MO.0.B.pcm",
    "TBX20_HUMAN.H11MO.0.D.pcm",
    "MEIS1_HUMAN.H11MO.1.B.pcm",
    "TGIF1_HUMAN.H11MO.0.A.pcm",
    "GATA4_HUMAN.H11MO.0.A.pcm",
    "TAL1_HUMAN.H11MO.0.A.pcm",
    "MAF_HUMAN.H11MO.1.B.pcm",
]
DEFAULT_PWM_NAMES_HUMAN = ["NKX25", "TBX20", "MEIS1", "TGIF1", "GATA4", "TAL1", "MAF"]

DEFAULT_PWM_MOUSE = [
    "NKX25_MOUSE.H11MO.0.A.pcm",
    "TBX20_MOUSE.H11MO.0.C.pcm",
    "MEIS1_MOUSE.H11MO.1.A.pcm",
    "TGIF1_MOUSE.H11MO.0.A.pcm",
    "HAND1_MOUSE.H11MO.0.C.pcm",
    "MAF_MOUSE.H11MO.1.A.pcm",
    "GATA1_MOUSE.H11MO.1.A.pcm",
    "GATA6_MOUSE.H11MO.0.A.pcm",
    "GATA4_MOUSE.H11MO.0.A.pcm",
]
DEFAULT_PWM_NAMES_MOUSE = ["Nkx25", "Tbx20", "Meis1", "Tgif1", "Hand1", "Maf", "Gata1", "Gata6", "Gata4"]


# ===========================================================================
# Argument parser construction
# ===========================================================================


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments shared by all subcommands."""
    parser.add_argument(
        "--organism", choices=["human", "mouse"], default="human",
        help="Organism. Determines default genome version and PWM files. (default: human)",
    )
    parser.add_argument(
        "--ensembl-release", type=int, default=None,
        help="Ensembl release number. Overrides organism default (75 for human/hg19, 54 for mouse/mm9).",
    )
    parser.add_argument(
        "--genome-name", type=str, default=None,
        help="Short genome identifier used for FASTA filenames (e.g. 'hg19', 'mm9'). "
             "Overrides organism default.",
    )
    parser.add_argument(
        "--path-fasta", default="",
        help="Directory containing per-chromosome FASTA files (e.g. hg19_chr1.fasta).",
    )
    parser.add_argument(
        "--genome-cache-dir", default=None,
        help="Directory for caching downloaded genome FASTAs. "
            "Default: ~/.cache/chipseq_pipeline/. "
            "Ignored when --path-fasta is set (legacy per-chromosome mode).",
    )
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity. (default: INFO)",
    )
    parser.add_argument(
        "--log-file", default=None,
        help="Write log output to this file in addition to stderr.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_pipeline",
        description="ChIP-seq peak analysis and TF binding-site detection pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    # ------------------------------------------------------------------
    # Subcommand: generador
    # ------------------------------------------------------------------
    p_gen = subparsers.add_parser(
        "generador",
        help="Load BED + RNA-seq data and generate binding-site and gene tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_common_args(p_gen)

    p_gen.add_argument("--bed-file", required=True,
                       help="BED filename without extension (located in --path-bed).")
    p_gen.add_argument("--rnaseq-file", required=True,
                       help="RNA-seq CSV filename without extension (located in --path-bed).")
    p_gen.add_argument("--output-prefix", required=True,
                       help="Base name for output files.")
    p_gen.add_argument("--path-bed", default="",
                       help="Directory containing BED and RNA-seq input files.")
    p_gen.add_argument("--path-out", default="",
                       help="Directory where output CSVs are written.")
    p_gen.add_argument("--path-pwm", default="",
                       help="Directory containing PWM/PCM files.")
    p_gen.add_argument("--pssm-file", default=None,
                       help="Primary TF PSSM filename (in --path-pwm). "
                            "Defaults to NKX25_HUMAN.H11MO.0.B.pcm (human) or "
                            "NKX25_MOUSE.H11MO.0.A.pcm (mouse).")
    p_gen.add_argument("--score-mult", type=float, default=0.9,
                       help="Fraction of PSSM max score used as the binding-site cutoff.")
    p_gen.add_argument("--binding-sites", nargs="*", default=None,
                       help="Confirmed binding-site sequences for list-based search. "
                            "Defaults to the 10 NKX2-5 core sequences.")
    p_gen.add_argument("--dist-max", type=int, default=100_000,
                       help="Maximum distance (bp) from peak to search for nearby genes.")
    p_gen.add_argument("--rnaseq-sep", default=";",
                       help="Field separator in the RNA-seq CSV.")
    p_gen.add_argument("--id-col", type=int, default=0,
                       help="Column index for gene ID in the RNA-seq file.")
    p_gen.add_argument("--updown-col", type=int, default=3,
                       help="Column index for log(fold-change) in the RNA-seq file.")
    p_gen.add_argument("--translate-col", type=int, default=None,
                       help="If set, this column contains gene names to translate to Ensembl IDs.")
    p_gen.add_argument("--test-mode", type=int, default=0,
                       help="If > 0, run on a random subset of this many peaks (for debugging).")

    # ------------------------------------------------------------------
    # Subcommand: memechip
    # ------------------------------------------------------------------
    p_meme = subparsers.add_parser(
        "memechip",
        help="Generate a FASTA file of peak/site sequences for MEME-ChIP.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_common_args(p_meme)

    p_meme.add_argument("--sites-file", required=True,
                        help="Input CSV filename without extension (peaks or binding sites).")
    p_meme.add_argument("--output-name", required=True,
                        help="Output FASTA filename without extension.")
    p_meme.add_argument("--path-sites", default="",
                        help="Directory containing the input CSV.")
    p_meme.add_argument("--path-out", default="",
                        help="Directory where the output FASTA is written.")
    p_meme.add_argument("--site-length", type=int, default=0,
                        help="If > 0, resize all sequences to this length (centred). "
                             "0 = keep original length.")
    p_meme.add_argument("--col-sites", nargs=3, type=int, default=[0, 1, 2],
                        metavar=("CHR_COL", "START_COL", "END_COL"),
                        help="Column indices for chr_n, pos_ini, pos_end in the input CSV.")
    p_meme.add_argument("--filter-col", type=int, default=None,
                        help="Column index to apply a value filter on.")
    p_meme.add_argument("--filter-val", default=None,
                        help="Value that triggers the filter (used with --filter-col).")
    p_meme.add_argument("--default-pass", action="store_true", default=True,
                        help="All sites pass by default; filter rule *excludes* matches. "
                             "Use --no-default-pass to invert (filter rule *includes* matches).")
    p_meme.add_argument("--no-default-pass", dest="default_pass", action="store_false")

    # ------------------------------------------------------------------
    # Subcommand: otros-tf
    # ------------------------------------------------------------------
    p_otf = subparsers.add_parser(
        "otros-tf",
        help="Find binding sites of additional TFs near primary-TF sites.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_common_args(p_otf)

    p_otf.add_argument("--sites-file", required=True,
                       help="Primary-TF sites CSV filename without extension.")
    p_otf.add_argument("--genes-file", required=True,
                       help="Nearby-genes CSV filename without extension.")
    p_otf.add_argument("--output-prefix", required=True,
                       help="Base name for output files.")
    p_otf.add_argument("--path-sites", default="",
                       help="Directory containing the primary-TF sites CSV.")
    p_otf.add_argument("--path-genes", default="",
                       help="Directory containing the genes CSV.")
    p_otf.add_argument("--path-out", default="",
                       help="Directory where output CSVs are written.")
    p_otf.add_argument("--path-pwm", default="",
                       help="Directory containing PWM/PCM files.")
    p_otf.add_argument("--pwm-files", nargs="+", default=None,
                       help="PWM filenames for additional TFs. "
                            "Defaults to the standard cardiac TF panel for the chosen organism.")
    p_otf.add_argument("--pwm-names", nargs="+", default=None,
                       help="Display names for the additional TFs (same order as --pwm-files). "
                            "Defaults to the standard panel names.")
    p_otf.add_argument("--dist-sites", type=int, default=250,
                       help="Distance (bp) to expand each primary-TF site when searching. "
                            "Use 0 to search only within the site.")
    p_otf.add_argument("--dist-max-genes", type=int, default=100_000,
                       help="Maximum distance to a gene +1 when reporting nearby genes.")
    p_otf.add_argument("--score-mult", type=float, default=0.9,
                       help="Fraction of PSSM max score used as the binding-site cutoff.")
    p_otf.add_argument("--score-report", action="store_true",
                       help="Print max-score and cutoff info for each PSSM before running "
                            "(useful for cutoff calibration).")

    return parser


# ===========================================================================
# Subcommand handlers
# ===========================================================================


def _resolve_genome(args) -> tuple:
    """Load the genome and return (genome, genome_name), applying any CLI overrides."""
    genome, genome_name = pipeline.load_genome(args.organism, ensembl_release=args.ensembl_release)
    if args.genome_name:
        genome_name = args.genome_name
    return genome, genome_name


def run_generador(args: argparse.Namespace) -> None:
    genome, genome_name = _resolve_genome(args)

    pssm_file = args.pssm_file or DEFAULT_PSSM[args.organism]
    pssm = pipeline.load_pssm(pssm_file, path_dir=args.path_pwm)
    score_cutoff = pssm.max * args.score_mult

    binding_sites = args.binding_sites if args.binding_sites is not None else DEFAULT_BINDING_SITES

    logging.info(
        "Running generador: organism=%s, genome=%s, pssm=%s, score_cutoff=%.4f, "
        "dist_max=%d, n_binding_sites=%d, test_mode=%d",
        args.organism, genome_name, pssm_file, score_cutoff,
        args.dist_max, len(binding_sites), args.test_mode,
    )

    pipeline.pipeline_generador(
        bed_file=args.bed_file,
        rnaseq_file=args.rnaseq_file,
        output_prefix=args.output_prefix,
        genome=genome,
        genome_name=genome_name,
        pssm=pssm,
        score_cutoff=score_cutoff,
        genome_cache_dir=args.genome_cache_dir,
        path_bed=args.path_bed,
        path_out=args.path_out,
        path_fasta=args.path_fasta,
        dist_max=args.dist_max,
        binding_sites=binding_sites,
        test_mode=args.test_mode,
        id_col_rnaseq=args.id_col,
        updown_col_rnaseq=args.updown_col,
        rnaseq_sep=args.rnaseq_sep,
        translate_col=args.translate_col,
    )


def run_memechip(args: argparse.Namespace) -> None:
    _, genome_name = _resolve_genome(args)

    filter_rules = []
    if args.filter_col is not None and args.filter_val is not None:
        filter_rules = [(args.filter_col, args.filter_val)]

    logging.info(
        "Running memechip: sites=%s, genome=%s, site_length=%d, n_filter_rules=%d",
        args.sites_file, genome_name, args.site_length, len(filter_rules),
    )

    pipeline.pipeline_meme_chip(
        sites_file=args.sites_file,
        output_name=args.output_name,
        genome_name=genome_name,
        site_length=args.site_length,
        genome_cache_dir=args.genome_cache_dir,
        col_sites=args.col_sites,
        filter_rules=filter_rules,
        default_pass_filter=args.default_pass,
        path_sites=args.path_sites,
        path_out=args.path_out,
        path_fasta=args.path_fasta,
    )


def run_otros_tf(args: argparse.Namespace) -> None:
    _, genome_name = _resolve_genome(args)

    # Default PWM panel per organism
    if args.pwm_files is None:
        pwm_files = DEFAULT_PWM_HUMAN if args.organism == "human" else DEFAULT_PWM_MOUSE
        pwm_names = DEFAULT_PWM_NAMES_HUMAN if args.organism == "human" else DEFAULT_PWM_NAMES_MOUSE
    else:
        pwm_files = args.pwm_files
        pwm_names = args.pwm_names or pwm_files  # fall back to filenames if names not provided

    logging.info(
        "Running otros-tf: sites=%s, genome=%s, n_pwm=%d, dist_sites=%d, score_mult=%.2f",
        args.sites_file, genome_name, len(pwm_files), args.dist_sites, args.score_mult,
    )

    pipeline.pipeline_otros_tf(
        sites_file=args.sites_file,
        genes_file=args.genes_file,
        output_prefix=args.output_prefix,
        genome_name=genome_name,
        pwm_files=pwm_files,
        dist_sites=args.dist_sites,
        dist_max_genes=args.dist_max_genes,
        genome_cache_dir=args.genome_cache_dir,
        pwm_names=pwm_names,
        path_sites=args.path_sites,
        path_genes=args.path_genes,
        path_out=args.path_out,
        path_fasta=args.path_fasta,
        path_pwm=args.path_pwm,
        score_mult=args.score_mult,
        score_report=args.score_report,
    )


# ===========================================================================
# Logging setup
# ===========================================================================


def _configure_logging(level: str, log_file: str | None) -> None:
    handlers = [logging.StreamHandler(sys.stderr)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
    )


# ===========================================================================
# Entry point
# ===========================================================================

HANDLERS = {
    "generador": run_generador,
    "memechip":  run_memechip,
    "otros-tf":  run_otros_tf,
}


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    _configure_logging(args.log_level, args.log_file)

    handler = HANDLERS.get(args.subcommand)
    if handler is None:
        parser.print_help()
        sys.exit(1)

    try:
        handler(args)
    except FileNotFoundError as exc:
        logging.error("File not found: %s", exc)
        sys.exit(1)
    except ValueError as exc:
        logging.error("Configuration error: %s", exc)
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("Interrupted by user.")
        sys.exit(0)


if __name__ == "__main__":
    main()