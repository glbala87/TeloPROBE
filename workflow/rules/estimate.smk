"""Rule: Estimator - Per-read, per-arm, per-sample TL summaries.

Input dependencies adapt based on skip_mapping and skip_clustering flags."""


def get_estimate_input(wildcards):
    """Determine estimate inputs based on pipeline mode."""
    base = f"{wildcards.output_dir}/{wildcards.sample}"
    inputs = {"segmented": f"{base}/segment/segmented.tsv"}

    if not config.get("skip_mapping", False):
        inputs["arm_assignments"] = f"{base}/anchor/arm_assignments.tsv"
        if not config.get("skip_clustering", False):
            inputs["clustered"] = f"{base}/cluster/allele_assignments.tsv"

    return inputs


def get_estimate_outputs(skip_mapping):
    """Return the correct set of output files."""
    outputs = {
        "read_level": "{output_dir}/{sample}/read_level.tsv",
        "sample_level": "{output_dir}/{sample}/sample_level.tsv",
    }
    # Always produce arm_level.tsv; when mapping is skipped it will be empty
    outputs["arm_level"] = "{output_dir}/{sample}/arm_level.tsv"
    return outputs


rule estimate:
    input:
        unpack(get_estimate_input)
    output:
        **get_estimate_outputs(config.get("skip_mapping", False))
    params:
        platform=config.get("platform", "ont"),
        mode=config.get("mode", "teloprobe"),
        bootstrap_n=config.get("bootstrap_n", 1000),
        bootstrap_alpha=config.get("bootstrap_alpha", 0.05),
        skip_mapping=config.get("skip_mapping", False),
    log:
        "{output_dir}/{sample}/logs/estimate.log"
    threads: 1
    run:
        import sys
        sys.path.insert(0, str(Path(workflow.basedir).parent / "src"))
        from teloprobe.estimate.arm_level import compute_arm_level
        from teloprobe.estimate.sample_level import compute_sample_level
        from teloprobe.constants import Platform, Mode
        import pandas as pd

        # Use the most complete data available
        if "clustered" in input.keys():
            read_level = pd.read_csv(input.clustered, sep="\t")
        elif "arm_assignments" in input.keys():
            segmented = pd.read_csv(input.segmented, sep="\t")
            arms = pd.read_csv(input.arm_assignments, sep="\t")
            if not arms.empty and "chr_arm" in arms.columns:
                arm_map = arms.set_index("read_id")["chr_arm"].to_dict()
                hap_map = arms.set_index("read_id")["haplotype"].to_dict()
                segmented["chromosome_arm"] = segmented["read_id"].map(arm_map).fillna("")
                segmented["haplotype"] = segmented["read_id"].map(hap_map).fillna("")
            read_level = segmented
        else:
            read_level = pd.read_csv(input.segmented, sep="\t")

        # Arm-level aggregation
        if params.skip_mapping:
            # No arm info available, write empty arm_level
            arm_level = pd.DataFrame(columns=[
                "chromosome_arm", "n_reads", "mean_tl", "median_tl",
                "ci_lower", "ci_upper"
            ])
        else:
            arm_level = compute_arm_level(
                read_level,
                bootstrap_n=params.bootstrap_n,
                bootstrap_alpha=params.bootstrap_alpha,
            )

        # Sample-level aggregation
        sample_level = compute_sample_level(
            read_level, arm_level,
            sample_id=wildcards.sample,
            platform=Platform(params.platform),
            mode=Mode(params.mode),
            bootstrap_n=params.bootstrap_n,
        )

        # Write outputs
        read_level.to_csv(output.read_level, sep="\t", index=False)
        arm_level.to_csv(output.arm_level, sep="\t", index=False)
        sample_level.to_csv(output.sample_level, sep="\t", index=False)
