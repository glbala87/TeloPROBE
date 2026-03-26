"""Rule: Anchor mapper - Subtelomeric alignment for arm assignment.

Reads segmented records from the segment pickle instead of re-reading
the original input file. This rule is only included when skip_mapping=False."""


rule anchor:
    input:
        segmented="{output_dir}/{sample}/segment/segmented.tsv",
        segmented_records="{output_dir}/{sample}/segment/segmented_records.pkl",
    output:
        arm_assignments="{output_dir}/{sample}/anchor/arm_assignments.tsv"
    params:
        reference=config.get("reference", ""),
        platform=config.get("platform", "ont"),
        identity_threshold=config.get("identity_threshold", 0.8),
        mapq_threshold=config.get("mapq_threshold", 20),
    log:
        "{output_dir}/{sample}/logs/anchor.log"
    threads:
        config.get("threads", 4)
    run:
        import sys, pickle
        sys.path.insert(0, str(Path(workflow.basedir).parent / "src"))
        from teloprobe.anchor.aligner import extract_flank_fastq, align_subtelomeric
        from teloprobe.anchor.arm_assigner import process_alignments
        from teloprobe.constants import Platform
        import pandas as pd

        segmented = pd.read_csv(input.segmented, sep="\t")
        ref = params.reference

        empty_df = pd.DataFrame(columns=[
            "read_id", "chr_arm", "haplotype",
            "alignment_identity", "mapping_quality", "anchor_confidence"
        ])

        if not ref:
            # No reference, output empty assignments
            empty_df.to_csv(output.arm_assignments, sep="\t", index=False)
        else:
            # Load segmented records from pickle
            with open(input.segmented_records, "rb") as fh:
                good_records = pickle.load(fh)

            sample_dir = Path(output.arm_assignments).parent
            good_ids = set(segmented[segmented["qc_status"] == "Good"]["read_id"]) if "qc_status" in segmented.columns else set(segmented["read_id"])
            bnd_map = dict(zip(
                segmented["read_id"],
                segmented["telomere_length_bp"]
            ))

            # Write flanks from serialized records
            flank_fq = str(sample_dir / "flanks.fastq")
            count = 0
            with open(flank_fq, "w") as fq:
                for rec in good_records:
                    if rec.read_id not in good_ids:
                        continue
                    bnd = bnd_map.get(rec.read_id)
                    if bnd is None or bnd <= 0:
                        continue
                    bnd = int(bnd)
                    flank = rec.sequence[bnd:]
                    if len(flank) < 200:
                        continue
                    fq.write(f"@{rec.read_id}\n{flank}\n+\n{'I'*len(flank)}\n")
                    count += 1

            if count > 0:
                aligned_bam = str(sample_dir / "aligned.bam")
                align_subtelomeric(
                    flank_fq, ref, aligned_bam,
                    threads=threads, platform=Platform(params.platform),
                )
                arm_df = process_alignments(
                    aligned_bam,
                    identity_threshold=params.identity_threshold,
                    mapq_threshold=params.mapq_threshold,
                )
                arm_df.to_csv(output.arm_assignments, sep="\t", index=False)
            else:
                empty_df.to_csv(output.arm_assignments, sep="\t", index=False)
