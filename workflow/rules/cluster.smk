"""Rule: Clusterer - Allele-specific clustering.

This rule is only included when both skip_mapping=False and
skip_clustering=False (controlled by conditional include in Snakefile)."""


rule cluster:
    input:
        segmented="{output_dir}/{sample}/segment/segmented.tsv",
        arm_assignments="{output_dir}/{sample}/anchor/arm_assignments.tsv"
    output:
        clustered="{output_dir}/{sample}/cluster/allele_assignments.tsv"
    params:
        allele_method=config.get("allele_method", "gmm"),
        allele_min_reads=config.get("allele_min_reads", 10),
    log:
        "{output_dir}/{sample}/logs/cluster.log"
    threads: 1
    run:
        import sys
        sys.path.insert(0, str(Path(workflow.basedir).parent / "src"))
        from teloprobe.cluster.allele_cluster import cluster_arms
        import pandas as pd

        segmented = pd.read_csv(input.segmented, sep="\t")
        arms = pd.read_csv(input.arm_assignments, sep="\t")

        # Merge arm assignments
        if not arms.empty and "chr_arm" in arms.columns:
            arm_map = arms.set_index("read_id")["chr_arm"].to_dict()
            hap_map = arms.set_index("read_id")["haplotype"].to_dict()
            segmented["chromosome_arm"] = segmented["read_id"].map(arm_map).fillna("")
            segmented["haplotype"] = segmented["read_id"].map(hap_map).fillna("")

        # Cluster
        result = cluster_arms(
            segmented,
            method=params.allele_method,
            min_reads=params.allele_min_reads,
        )

        result.to_csv(output.clustered, sep="\t", index=False)
