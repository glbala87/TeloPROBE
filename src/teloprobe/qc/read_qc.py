"""Per-read QC summary and diagnostics."""

import pandas as pd
import numpy as np


def qc_summary_table(read_level: pd.DataFrame) -> pd.DataFrame:
    """Generate QC summary table with counts and statistics per status.

    Returns DataFrame with columns:
    - status: QC status label
    - total_reads: count of reads with this status
    - median_read_length: median read length
    - median_tl: median telomere length (for Good reads)
    - fraction: fraction of total reads
    """
    rows = []
    total = len(read_level)

    for status, group in read_level.groupby("qc_status"):
        n = len(group)
        tl_vals = group["telomere_length_bp"].values
        tl_positive = tl_vals[tl_vals > 0] if len(tl_vals) > 0 else np.array([])

        rows.append({
            "status": status,
            "total_reads": n,
            "median_read_length": int(group["read_length"].median()),
            "median_tl": float(np.median(tl_positive)) if len(tl_positive) > 0 else 0.0,
            "fraction": round(n / total, 4) if total > 0 else 0.0,
        })

    df = pd.DataFrame(rows).sort_values("total_reads", ascending=False)
    return df.reset_index(drop=True)


def filter_breakdown(read_level: pd.DataFrame) -> dict:
    """Return a dict of QC status -> count for waterfall/funnel charts."""
    return read_level["qc_status"].value_counts().to_dict()


def terminality_failure_rate(read_level: pd.DataFrame) -> float:
    """Fraction of reads that failed due to non-terminal telomere position."""
    terminal_failures = {"TooCloseStart", "TooCloseEnd", "StartNotRepeats"}
    total = len(read_level)
    if total == 0:
        return 0.0
    failed = read_level["qc_status"].isin(terminal_failures).sum()
    return float(failed / total)


def internal_telomere_rate(read_level: pd.DataFrame) -> float:
    """Fraction of reads with telomere-side classified as 'internal'."""
    if "telomere_side" not in read_level.columns:
        return 0.0
    total = len(read_level)
    if total == 0:
        return 0.0
    internal = (read_level["telomere_side"] == "internal").sum()
    return float(internal / total)
