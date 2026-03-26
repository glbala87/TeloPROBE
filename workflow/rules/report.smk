"""Rule: Report - HTML report generation.

Fixed to properly declare all per-sample inputs and build a single
multi-sample report."""


def get_report_inputs(wildcards):
    """Gather all sample outputs for multi-sample report.

    Explicitly declares every file the report rule reads so that
    Snakemake tracks all dependencies correctly.
    """
    inputs = {}
    for i, sample in enumerate(SAMPLES):
        base = f"{wildcards.output_dir}/{sample}"
        inputs[f"read_{i}"] = f"{base}/read_level.tsv"
        inputs[f"sample_{i}"] = f"{base}/sample_level.tsv"
        inputs[f"arm_{i}"] = f"{base}/arm_level.tsv"
        inputs[f"qc_{i}"] = f"{base}/qc_summary.tsv"
    return inputs


rule report:
    input:
        unpack(get_report_inputs)
    output:
        report="{output_dir}/report.html"
    params:
        platform=config.get("platform", "ont"),
        mode=config.get("mode", "teloprobe"),
        samples=SAMPLES,
    log:
        "{output_dir}/logs/report.log"
    threads: 1
    run:
        import sys
        sys.path.insert(0, str(Path(workflow.basedir).parent / "src"))
        from teloprobe.config import Config
        from teloprobe.report.builder import build_report
        from teloprobe.constants import Platform, Mode
        import pandas as pd

        # Collect data across all samples
        all_read_level = []
        all_arm_level = []
        all_sample_level = []

        for sample in params.samples:
            base = Path(f"{wildcards.output_dir}/{sample}")

            rl = pd.read_csv(base / "read_level.tsv", sep="\t")
            rl["sample_id"] = sample
            all_read_level.append(rl)

            al = pd.read_csv(base / "arm_level.tsv", sep="\t")
            if not al.empty:
                al["sample_id"] = sample
            all_arm_level.append(al)

            sl = pd.read_csv(base / "sample_level.tsv", sep="\t")
            all_sample_level.append(sl)

        read_level = pd.concat(all_read_level, ignore_index=True)
        arm_level = pd.concat(all_arm_level, ignore_index=True) if all_arm_level else pd.DataFrame()
        sample_level = pd.concat(all_sample_level, ignore_index=True) if all_sample_level else pd.DataFrame()

        cfg = Config()
        cfg.platform = Platform(params.platform)
        cfg.mode = Mode(params.mode)
        # For single sample, use the sample name; for multi, use "multi"
        if len(params.samples) == 1:
            cfg.sample_id = params.samples[0]
        else:
            cfg.sample_id = "multi"

        build_report(
            read_level, arm_level, sample_level,
            output.report, cfg,
        )
