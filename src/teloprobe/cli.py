"""TeloPROBE CLI entry point using Click."""

import logging
import sys
import time
from pathlib import Path

import click
import numpy as np

from .constants import Mode, Platform
from .config import Config


def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


@click.group()
@click.version_option(package_name="teloprobe")
def cli():
    """TeloPROBE: Robust telomere length estimation from long-read sequencing."""
    pass


@cli.command()
@click.option("-i", "--input", "input_path", required=True,
              help="Input BAM/CRAM/FASTQ file or directory.")
@click.option("-o", "--output", "output_dir", default="output",
              help="Output directory.")
@click.option("-r", "--reference", default=None,
              help="Subtelomeric reference FASTA (for arm assignment).")
@click.option("--mode", type=click.Choice(["wgs", "teloprobe"]), default="teloprobe",
              help="Analysis mode: wgs (whole-genome) or teloprobe (enriched).")
@click.option("--platform", type=click.Choice(["ont", "pacbio", "auto"]),
              default="auto", help="Sequencing platform.")
@click.option("--sample-id", default="sample",
              help="Sample identifier.")
@click.option("--sample-sheet", default=None,
              help="CSV sample sheet for multi-sample runs.")
@click.option("-t", "--threads", default=4, type=int,
              help="Number of threads.")
@click.option("--skip-mapping", is_flag=True,
              help="Skip alignment and arm assignment.")
@click.option("--skip-clustering", is_flag=True,
              help="Skip allele clustering.")
@click.option("--config", "config_path", default=None,
              help="Custom YAML configuration file.")
@click.option("--min-repeats", default=100, type=int,
              help="Minimum telomeric repeat count.")
@click.option("--min-length", default=100, type=int,
              help="Minimum read length.")
@click.option("--min-quality", default=9, type=float,
              help="Minimum mean read quality.")
@click.option("--changepoint-method", type=click.Choice(["pelt", "binseg", "bayesian"]),
              default="pelt", help="Change-point detection method.")
@click.option("--bootstrap-n", default=1000, type=int,
              help="Number of bootstrap resamples for CIs.")
@click.option("--offline", is_flag=True,
              help="Embed Plotly JS inline for offline HTML reports (requires plotly).")
@click.option("-v", "--verbose", is_flag=True, help="Verbose output.")
def run(input_path, output_dir, reference, mode, platform, sample_id,
        sample_sheet, threads, skip_mapping, skip_clustering, config_path,
        min_repeats, min_length, min_quality, changepoint_method,
        bootstrap_n, offline, verbose):
    """Run the full TeloPROBE analysis pipeline."""
    setup_logging(verbose)
    logger = logging.getLogger("teloprobe")

    pipeline_start = time.time()

    from .config import Config

    # Load or create config
    if config_path:
        config = Config.from_yaml(config_path)
    else:
        config = Config()

    # Apply CLI overrides
    config.output_dir = output_dir
    config.sample_id = sample_id
    config.mode = Mode(mode)
    config.threads = threads
    config.skip_mapping = skip_mapping or reference is None
    config.min_repeats = min_repeats
    config.min_length = min_length
    config.min_quality = int(min_quality)
    config.changepoint_method = changepoint_method
    config.bootstrap_n = bootstrap_n
    config.offline = offline

    if offline:
        try:
            import plotly.offline  # noqa: F401
        except ImportError:
            click.echo(
                "WARNING: --offline requires the 'plotly' package. "
                "Install it with: pip install plotly\n"
                "Falling back to CDN-based Plotly.",
                err=True,
            )

    if reference:
        config.reference = reference

    # ------------------------------------------------------------------
    # Sample-sheet mode: loop over each sample in the sheet
    # ------------------------------------------------------------------
    if sample_sheet:
        from .ingest.sample_sheet import parse_sample_sheet, validate_sample_sheet

        logger.info("Parsing sample sheet: %s", sample_sheet)
        samples = parse_sample_sheet(sample_sheet)

        errors = validate_sample_sheet(samples)
        if errors:
            for err in errors:
                logger.error("Sample sheet error: %s", err)
            sys.exit(1)

        logger.info("Found %d sample(s) in sample sheet", len(samples))

        import pandas as pd

        all_sample_dfs = []

        for idx, sc in enumerate(samples, start=1):
            sample_start = time.time()
            logger.info(
                "=== Processing sample %d/%d: %s ===",
                idx, len(samples), sc.sample_id,
            )

            # Build a per-sample config
            sample_config = Config.from_yaml(config_path) if config_path else Config()
            sample_config.input_path = sc.input_path
            sample_config.output_dir = output_dir
            sample_config.sample_id = sc.sample_id
            sample_config.mode = sc.mode
            sample_config.platform = sc.platform
            sample_config.threads = threads
            sample_config.skip_mapping = skip_mapping or (
                sc.reference is None and reference is None
            )
            sample_config.min_repeats = min_repeats
            sample_config.min_length = min_length
            sample_config.min_quality = int(min_quality)
            sample_config.changepoint_method = changepoint_method
            sample_config.bootstrap_n = bootstrap_n

            if sc.reference:
                sample_config.reference = sc.reference
            elif reference:
                sample_config.reference = reference

            sample_config.apply_platform_preset()

            sample_level_df = _run_pipeline(sample_config, skip_clustering)
            if sample_level_df is not None and not sample_level_df.empty:
                all_sample_dfs.append(sample_level_df)

            elapsed = time.time() - sample_start
            logger.info(
                "Sample %s completed in %.1f seconds", sc.sample_id, elapsed,
            )

        # Aggregate all sample-level results into a combined file
        if all_sample_dfs:
            combined = pd.concat(all_sample_dfs, ignore_index=True)
            out = Path(output_dir)
            out.mkdir(parents=True, exist_ok=True)
            combined.to_csv(out / "all_samples.tsv", sep="\t", index=False)
            logger.info(
                "Aggregated results for %d samples written to %s",
                len(all_sample_dfs), out / "all_samples.tsv",
            )

        total_elapsed = time.time() - pipeline_start
        logger.info("All samples completed in %.1f seconds", total_elapsed)
        return

    # ------------------------------------------------------------------
    # Single-sample mode
    # ------------------------------------------------------------------
    config.input_path = input_path

    # Auto-detect platform
    if platform == "auto":
        from .ingest.platform_detect import detect_platform
        config.platform = detect_platform(input_path)
        if config.platform == Platform.UNKNOWN:
            config.platform = Platform.ONT
            logger.info("Could not detect platform, defaulting to ONT")
    else:
        config.platform = Platform(platform)

    config.apply_platform_preset()

    logger.info("TeloPROBE v2.0.0 starting")
    logger.info("  Input: %s", input_path)
    logger.info("  Mode: %s | Platform: %s", config.mode.value, config.platform.value)
    logger.info("  Output: %s", output_dir)

    # Run pipeline
    _run_pipeline(config, skip_clustering)

    total_elapsed = time.time() - pipeline_start
    logger.info("Pipeline completed in %.1f seconds", total_elapsed)


def _run_pipeline(config: Config, skip_clustering: bool = False):
    """Execute the analysis pipeline for a single sample.

    Returns the sample-level DataFrame (useful for multi-sample aggregation).

    NOTE on parallelism: The per-read HMM/PELT computations are CPU-bound.
    In principle a ``concurrent.futures.ProcessPoolExecutor`` could speed up
    ``compute_read_level``, but hmmlearn model objects are not always
    pickle-safe across platforms, and forking with NumPy/OpenBLAS can
    deadlock.  For now, parallelism is limited to the alignment step
    (minimap2 uses ``config.threads`` internally).  If you need per-read
    parallelism, consider running separate processes per sample via the
    sample-sheet mode, which avoids pickling model objects altogether.
    """
    import pandas as pd
    from pathlib import Path

    logger = logging.getLogger("teloprobe")
    out = Path(config.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Module 1: Ingest
    # ------------------------------------------------------------------
    t0 = time.time()
    logger.info("Step 1/8: Ingesting reads...")
    from .ingest.reader import read_input
    records = list(read_input(
        config.input_path,
        platform=config.platform,
        mode=config.mode,
        sample_id=config.sample_id,
        min_length=config.min_length,
        min_quality=config.min_quality,
    ))
    logger.info(
        "  Loaded %d reads (%.1f s)", len(records), time.time() - t0,
    )

    if not records:
        logger.error("No reads loaded. Check input path and filters.")
        click.echo(
            f"\nTeloPROBE: 0 reads passed filters for sample "
            f"'{config.sample_id}'. Nothing to analyse."
        )
        # Write a minimal sample-level output so downstream tooling
        # can still parse the results directory.
        sample_dir = out / config.sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        empty_sample = pd.DataFrame([{
            "sample_id": config.sample_id,
            "total_reads": 0,
            "informative_reads": 0,
            "median_tl": float("nan"),
            "global_ci_lower": float("nan"),
            "global_ci_upper": float("nan"),
            "arms_assigned": 0,
            "qc_status": "FAIL_NO_READS",
        }])
        empty_sample.to_csv(
            sample_dir / "sample_level.tsv", sep="\t", index=False,
        )
        return empty_sample

    # ------------------------------------------------------------------
    # Modules 2-3: Candidate finding + Segmentation
    # ------------------------------------------------------------------
    t0 = time.time()
    logger.info("Step 2/8: Scanning motifs and detecting boundaries...")
    from .estimate.read_level import compute_read_level
    read_level_df = compute_read_level(records, config)
    logger.info("  Read-level analysis done (%.1f s)", time.time() - t0)

    # Check whether any reads survived the analysis steps
    informative = read_level_df[read_level_df["qc_status"] == "Good"]
    if informative.empty:
        logger.warning(
            "0 reads passed QC for sample '%s'. "
            "Writing empty outputs and skipping remaining steps.",
            config.sample_id,
        )
        sample_dir = out / config.sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        from .io.writers import write_read_level
        write_read_level(read_level_df, sample_dir / "read_level.tsv")
        empty_sample = pd.DataFrame([{
            "sample_id": config.sample_id,
            "total_reads": len(read_level_df),
            "informative_reads": 0,
            "median_tl": float("nan"),
            "global_ci_lower": float("nan"),
            "global_ci_upper": float("nan"),
            "arms_assigned": 0,
            "qc_status": "FAIL_NO_INFORMATIVE_READS",
        }])
        empty_sample.to_csv(
            sample_dir / "sample_level.tsv", sep="\t", index=False,
        )
        click.echo(
            f"\nTeloPROBE: 0 informative reads for sample "
            f"'{config.sample_id}'. Outputs written but no estimates made."
        )
        return empty_sample

    # ------------------------------------------------------------------
    # Module 4: Anchor mapping (optional)
    # ------------------------------------------------------------------
    arm_level_df = pd.DataFrame()
    if not config.skip_mapping and config.reference:
        t0 = time.time()
        logger.info("Step 4/8: Aligning subtelomeric flanks...")
        arm_level_df = _run_alignment(read_level_df, records, config, out)
        logger.info("  Alignment done (%.1f s)", time.time() - t0)
    else:
        logger.info("Step 4/8: Skipping alignment (no reference or --skip-mapping)")

    # ------------------------------------------------------------------
    # Module 5: Allele clustering (optional)
    # ------------------------------------------------------------------
    if not skip_clustering and not arm_level_df.empty:
        t0 = time.time()
        logger.info("Step 5/8: Clustering alleles...")
        from .cluster.allele_cluster import cluster_arms
        read_level_df = cluster_arms(
            read_level_df, method=config.allele_method,
            min_reads=config.allele_min_reads,
        )
        logger.info("  Clustering done (%.1f s)", time.time() - t0)

    # ------------------------------------------------------------------
    # Module 6: Estimation
    # ------------------------------------------------------------------
    t0 = time.time()
    logger.info("Step 6/8: Computing estimates...")
    from .estimate.arm_level import compute_arm_level
    from .estimate.sample_level import compute_sample_level

    # Check whether any reads have an arm assignment.  The column may
    # contain empty strings, NaN, or None for unassigned reads.
    has_arm = (
        "chromosome_arm" in read_level_df.columns
        and read_level_df["chromosome_arm"]
            .replace("", np.nan)
            .dropna()
            .shape[0] > 0
    )

    if has_arm:
        arm_level_df = compute_arm_level(
            read_level_df,
            bootstrap_n=config.bootstrap_n,
            bootstrap_alpha=config.bootstrap_alpha,
        )

    sample_level_df = compute_sample_level(
        read_level_df, arm_level_df,
        sample_id=config.sample_id,
        platform=config.platform,
        mode=config.mode,
        bootstrap_n=config.bootstrap_n,
    )
    logger.info("  Estimation done (%.1f s)", time.time() - t0)

    # ------------------------------------------------------------------
    # Module 7: QC
    # ------------------------------------------------------------------
    t0 = time.time()
    logger.info("Step 7/8: Running QC checks...")
    from .qc.sample_qc import sample_qc_checks, flag_low_coverage_arms
    from .qc.validators import cross_validate_evidence
    from .qc.read_qc import qc_summary_table
    from .qc.sample_qc import downsampling_stability

    qc_checks = sample_qc_checks(read_level_df, arm_level_df)
    read_level_df = cross_validate_evidence(read_level_df)

    if not arm_level_df.empty:
        arm_level_df = flag_low_coverage_arms(arm_level_df)

    # Downsampling stability
    good_tl = read_level_df[
        read_level_df["qc_status"] == "Good"
    ]["telomere_length_bp"].values.astype(float)
    good_tl = good_tl[good_tl > 0]
    stability = pd.DataFrame()
    if len(good_tl) >= 50:
        stability = downsampling_stability(good_tl, random_seed=config.random_seed)
    logger.info("  QC done (%.1f s)", time.time() - t0)

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    t0 = time.time()
    logger.info("Writing outputs...")
    from .io.writers import (
        write_read_level, write_arm_level, write_sample_level,
        write_qc_json, write_config_json,
    )

    sample_dir = out / config.sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    write_read_level(read_level_df, sample_dir / "read_level.tsv")
    write_arm_level(arm_level_df, sample_dir / "arm_level.tsv")
    write_sample_level(sample_level_df, sample_dir / "sample_level.tsv")
    write_qc_json(qc_checks, sample_dir / "qc_flags.json")
    write_config_json(config, out / "params.json")

    # ------------------------------------------------------------------
    # Module 8: Report
    # ------------------------------------------------------------------
    logger.info("Step 8/8: Generating report...")
    from .report.builder import build_report
    build_report(
        read_level_df, arm_level_df, sample_level_df,
        out / "report.html",
        config,
        qc_checks=qc_checks,
        stability_data=stability,
    )
    logger.info("  Output written (%.1f s)", time.time() - t0)

    # Print summary
    if not sample_level_df.empty:
        row = sample_level_df.iloc[0]
        click.echo(f"\nTeloPROBE Analysis Complete")
        click.echo(f"  Sample: {config.sample_id}")
        click.echo(f"  Total reads: {row.get('total_reads', 0)}")
        click.echo(f"  Informative reads: {row.get('informative_reads', 0)}")
        click.echo(f"  Median TL: {row.get('median_tl', 0):.0f} bp")
        click.echo(f"  95% CI: [{row.get('global_ci_lower', 0):.0f}, {row.get('global_ci_upper', 0):.0f}] bp")
        if row.get('arms_assigned', 0) > 0:
            click.echo(f"  Arms assigned: {row.get('arms_assigned', 0)}")
        click.echo(f"  Report: {out / 'report.html'}")

    return sample_level_df


def _run_alignment(read_level_df, records, config, out):
    """Run subtelomeric alignment and arm assignment."""
    import tempfile
    import pandas as pd

    logger = logging.getLogger("teloprobe")
    sample_dir = out / config.sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    # Write good reads to temporary FASTQ
    good_reads = read_level_df[read_level_df["qc_status"] == "Good"]
    if good_reads.empty:
        return pd.DataFrame()

    # Create boundaries TSV for flank extraction
    boundaries_path = sample_dir / "boundaries.tsv"
    boundaries = read_level_df[["read_id", "telomere_length_bp", "qc_status"]].copy()
    boundaries.columns = ["read_id", "boundary_pos", "qc_status"]
    boundaries.to_csv(boundaries_path, sep="\t", index=False)

    # Write reads to temporary BAM
    import pysam
    tmp_bam = sample_dir / "candidates.bam"

    # For now, create a FASTQ of subtelomeric flanks directly
    flank_fastq = sample_dir / "flanks.fastq"
    record_map = {r.read_id: r for r in records}
    count = 0

    with open(flank_fastq, "w") as fq:
        for _, row in good_reads.iterrows():
            read_id = row["read_id"]
            tl = int(row["telomere_length_bp"])
            if tl <= 0 or read_id not in record_map:
                continue

            rec = record_map[read_id]
            flank = rec.sequence[tl:]
            if len(flank) < 200:
                continue

            # Convert quality scores from the ReadRecord numpy array.
            # rec.qualities is a numpy array of Phred scores (integers).
            # We need to produce a FASTQ quality string (Phred+33 ASCII).
            if rec.qualities is not None and len(rec.qualities) > tl:
                flank_quals = rec.qualities[tl:]
                # Clip to valid ASCII range: Phred 0-93 maps to '!' (33)
                # through '~' (126).
                clipped = np.clip(flank_quals, 0, 93).astype(np.int32) + 33
                qual_str = "".join(chr(q) for q in clipped)
            else:
                # Fallback: use Phred 20 ('5') as a reasonable default
                # when quality information is genuinely unavailable.
                qual_str = "5" * len(flank)

            fq.write(f"@{read_id}\n{flank}\n+\n{qual_str}\n")
            count += 1

    if count == 0:
        logger.warning("No subtelomeric flanks to align")
        return pd.DataFrame()

    logger.info("  Extracted %d subtelomeric flanks", count)

    # Align
    from .anchor.aligner import align_subtelomeric
    aligned_bam = sample_dir / "aligned.bam"
    align_subtelomeric(
        flank_fastq, config.reference, aligned_bam,
        threads=config.threads, platform=config.platform,
    )

    # Process alignments
    from .anchor.arm_assigner import process_alignments
    arm_df = process_alignments(
        str(aligned_bam),
        identity_threshold=config.identity_threshold,
        mapq_threshold=config.mapq_threshold,
    )

    # Merge arm assignments back into read_level
    if not arm_df.empty:
        arm_map = arm_df.set_index("read_id")[
            ["chr_arm", "haplotype", "alignment_identity",
             "mapping_quality", "anchor_confidence"]
        ].to_dict("index")

        for idx, row in read_level_df.iterrows():
            rid = row["read_id"]
            if rid in arm_map:
                info = arm_map[rid]
                read_level_df.at[idx, "chromosome_arm"] = info.get("chr_arm", "")
                read_level_df.at[idx, "haplotype"] = info.get("haplotype", "")
                read_level_df.at[idx, "alignment_identity"] = info.get("alignment_identity", "")
                read_level_df.at[idx, "mapping_quality"] = info.get("mapping_quality", "")
                read_level_df.at[idx, "anchor_confidence"] = info.get("anchor_confidence", "none")

    return arm_df


@cli.command()
@click.argument("results_dir")
@click.option("-o", "--output", default="report.html",
              help="Output HTML report path.")
def report(results_dir, output):
    """Re-generate HTML report from existing results."""
    setup_logging()
    import pandas as pd
    from .config import Config
    from .report.builder import build_report

    results = Path(results_dir)

    read_level = pd.read_csv(results / "read_level.tsv", sep="\t")
    arm_level_path = results / "arm_level.tsv"
    arm_level = pd.read_csv(arm_level_path, sep="\t") if arm_level_path.exists() else pd.DataFrame()
    sample_level_path = results / "sample_level.tsv"
    sample_level = pd.read_csv(sample_level_path, sep="\t") if sample_level_path.exists() else pd.DataFrame()

    config = Config()
    build_report(read_level, arm_level, sample_level, output, config)
    click.echo(f"Report written to {output}")


@cli.command()
@click.option("-i", "--input", "input_path", required=True)
def validate(input_path):
    """Validate input files without running analysis."""
    setup_logging()
    from .ingest.reader import detect_format

    path = Path(input_path)
    if not path.exists():
        click.echo(f"ERROR: File not found: {path}")
        sys.exit(1)

    fmt = detect_format(path)
    click.echo(f"File: {path}")
    click.echo(f"Format: {fmt}")

    from .ingest.platform_detect import detect_platform
    if fmt in ("bam", "cram"):
        platform = detect_platform(path)
        click.echo(f"Platform: {platform.value}")

    click.echo("Validation passed.")


if __name__ == "__main__":
    cli()
