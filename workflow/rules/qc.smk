"""Rule: QC - Quality control and validation.

Adapts inputs based on whether arm_level data is available."""


rule qc:
    input:
        read_level="{output_dir}/{sample}/read_level.tsv",
        arm_level="{output_dir}/{sample}/arm_level.tsv",
    output:
        qc_summary="{output_dir}/{sample}/qc_summary.tsv",
        qc_flags="{output_dir}/{sample}/qc_flags.json"
    log:
        "{output_dir}/{sample}/logs/qc.log"
    threads: 1
    run:
        import sys
        sys.path.insert(0, str(Path(workflow.basedir).parent / "src"))
        from teloprobe.qc.read_qc import qc_summary_table
        from teloprobe.qc.sample_qc import sample_qc_checks
        from teloprobe.qc.validators import cross_validate_evidence
        from teloprobe.io.writers import write_qc_json
        import pandas as pd

        read_level = pd.read_csv(input.read_level, sep="\t")
        arm_level = pd.read_csv(input.arm_level, sep="\t")

        # QC summary table
        qc_summary = qc_summary_table(read_level)
        qc_summary.to_csv(output.qc_summary, sep="\t", index=False)

        # QC checks (handle empty arm_level gracefully)
        if arm_level.empty:
            checks = sample_qc_checks(read_level, pd.DataFrame())
        else:
            checks = sample_qc_checks(read_level, arm_level)
        write_qc_json(checks, output.qc_flags)

        # Cross-validation (updates read_level in place for report)
        cross_validate_evidence(read_level)
