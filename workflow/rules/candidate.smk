"""Rule: Candidate finder - Fast motif scan for telomeric reads.

Reads from the ingest pickle instead of re-reading the original input file."""


rule candidate_scan:
    input:
        reads="{output_dir}/{sample}/ingest/reads.tsv",
        records="{output_dir}/{sample}/ingest/records.pkl",
    output:
        candidates="{output_dir}/{sample}/candidate/candidates.tsv",
        candidate_records="{output_dir}/{sample}/candidate/candidate_records.pkl",
    params:
        min_repeats=config.get("min_repeats", 100),
        mode=config.get("mode", "teloprobe"),
    log:
        "{output_dir}/{sample}/logs/candidate.log"
    threads: 1
    run:
        import sys, pickle
        sys.path.insert(0, str(Path(workflow.basedir).parent / "src"))
        from teloprobe.candidate.motif_scanner import quick_filter, scan_motifs, telomeric_fraction
        from teloprobe.candidate.read_filter import trim_adapters
        from teloprobe.constants import Mode
        import pandas as pd

        # Load serialized records from ingest
        with open(input.records, "rb") as fh:
            all_records = pickle.load(fh)

        rows = []
        candidate_records = []
        for rec in all_records:
            if Mode(params.mode) == Mode.TELOSEQ:
                rec, _ = trim_adapters(rec)

            is_candidate = quick_filter(rec.sequence, min_repeats=params.min_repeats)

            if is_candidate:
                motif_array = scan_motifs(rec.sequence)
                rows.append({
                    "read_id": rec.read_id,
                    "read_length": rec.read_length,
                    "motif_density": telomeric_fraction(motif_array),
                    "is_candidate": True,
                })
                candidate_records.append(rec)

        df = pd.DataFrame(rows)
        df.to_csv(output.candidates, sep="\t", index=False)

        # Serialize candidate ReadRecords for downstream
        with open(output.candidate_records, "wb") as fh:
            pickle.dump(candidate_records, fh, protocol=pickle.HIGHEST_PROTOCOL)
