"""Rule: Ingest - Read BAM/CRAM/FASTQ inputs, filter, and serialize to
an intermediate pickle so downstream rules never re-read the original file."""


def get_input_path(wildcards):
    """Resolve input path for a sample."""
    sample_inputs = config.get("sample_inputs", {})
    if wildcards.sample in sample_inputs:
        return sample_inputs[wildcards.sample]
    return config.get("input_path", "")


rule ingest:
    input:
        get_input_path
    output:
        reads="{output_dir}/{sample}/ingest/reads.tsv",
        records="{output_dir}/{sample}/ingest/records.pkl",
    params:
        platform=config.get("platform", "ont"),
        mode=config.get("mode", "teloprobe"),
        min_length=config.get("min_length", 100),
        min_quality=config.get("min_quality", 9),
    log:
        "{output_dir}/{sample}/logs/ingest.log"
    threads: 1
    run:
        import sys, pickle
        sys.path.insert(0, str(Path(workflow.basedir).parent / "src"))
        from teloprobe.ingest.reader import read_input
        from teloprobe.constants import Platform, Mode
        import pandas as pd

        records = []
        rows = []
        for rec in read_input(
            input[0],
            platform=Platform(params.platform),
            mode=Mode(params.mode),
            sample_id=wildcards.sample,
            min_length=params.min_length,
            min_quality=params.min_quality,
        ):
            records.append(rec)
            rows.append({
                "read_id": rec.read_id,
                "read_length": rec.read_length,
                "mean_quality": rec.mean_quality,
            })

        # Write summary TSV (for quick inspection)
        df = pd.DataFrame(rows)
        df.to_csv(output.reads, sep="\t", index=False)

        # Serialize full ReadRecord objects for downstream rules
        with open(output.records, "wb") as fh:
            pickle.dump(records, fh, protocol=pickle.HIGHEST_PROTOCOL)
