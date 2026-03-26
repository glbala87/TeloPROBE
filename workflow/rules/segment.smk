"""Rule: Segmenter - HMM/change-point boundary detection.

Reads candidate records from the candidate pickle instead of re-reading
the original input file."""


rule segment:
    input:
        candidates="{output_dir}/{sample}/candidate/candidates.tsv",
        candidate_records="{output_dir}/{sample}/candidate/candidate_records.pkl",
    output:
        segmented="{output_dir}/{sample}/segment/segmented.tsv",
        segmented_records="{output_dir}/{sample}/segment/segmented_records.pkl",
    params:
        platform=config.get("platform", "ont"),
        mode=config.get("mode", "teloprobe"),
        changepoint_method=config.get("changepoint_method", "pelt"),
        random_seed=config.get("random_seed", 42),
    log:
        "{output_dir}/{sample}/logs/segment.log"
    threads: 1
    run:
        import sys, pickle
        sys.path.insert(0, str(Path(workflow.basedir).parent / "src"))
        from teloprobe.config import Config
        from teloprobe.constants import Platform, Mode
        from teloprobe.estimate.read_level import compute_read_level
        import pandas as pd

        # Load candidate records from candidate pickle
        with open(input.candidate_records, "rb") as fh:
            records = pickle.load(fh)

        # Build config
        cfg = Config()
        cfg.platform = Platform(params.platform)
        cfg.mode = Mode(params.mode)
        cfg.changepoint_method = params.changepoint_method
        cfg.random_seed = params.random_seed
        cfg.apply_platform_preset()

        # Process reads
        df = compute_read_level(records, cfg)
        df.to_csv(output.segmented, sep="\t", index=False)

        # Build a lookup of records that passed segmentation (Good QC)
        good_ids = set(df[df.get("qc_status", pd.Series(dtype=str)) == "Good"]["read_id"]) if "qc_status" in df.columns else set(df["read_id"])
        good_records = [r for r in records if r.read_id in good_ids]

        with open(output.segmented_records, "wb") as fh:
            pickle.dump(good_records, fh, protocol=pickle.HIGHEST_PROTOCOL)
