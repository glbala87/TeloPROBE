"""Wrapper for Telogator2 (Stephens & Kocher, 2024).

Telogator2 performs allele-specific telomere length estimation via TVR-based
clustering. It bundles its own subtelomeric reference and requires minimap2.

GitHub: https://github.com/zstephens/telogator2
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .base import ToolWrapper, ToolResult

logger = logging.getLogger(__name__)


class Telogator2Wrapper(ToolWrapper):
    """Wrapper for Telogator2."""

    def __init__(self):
        super().__init__("Telogator2")

    def check_installed(self) -> bool:
        """Telogator2 can be installed as CLI or run as script."""
        if self._check_binary("telogator2"):
            return True
        # Check for script-based installation
        rc, _, _, _ = self._run_command(
            ["python", "-c", "import telogator2"], timeout=10
        )
        return rc == 0

    def run(
        self,
        input_path: str,
        output_dir: str,
        platform: str = "ont",
        reference: Optional[str] = None,
        threads: int = 4,
        extra_args: Optional[dict] = None,
    ) -> ToolResult:
        """Run Telogator2.

        Args:
            platform: "ont" or "pacbio" → mapped to Telogator2's -r {ont, hifi}
            reference: Not needed — Telogator2 bundles its own reference.
        """
        extra_args = extra_args or {}
        read_type = "hifi" if platform == "pacbio" else "ont"

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        cmd = ["telogator2"]
        # Fall back to python script if CLI not on PATH
        if not self._check_binary("telogator2"):
            cmd = ["python", "-m", "telogator2"]

        cmd.extend([
            "-i", input_path,
            "-o", output_dir,
            "-r", read_type,
            "-p", str(threads),
        ])

        # Optional: minimum read length
        if "min_length" in extra_args:
            cmd.extend(["-l", str(extra_args["min_length"])])
        else:
            cmd.extend(["-l", "4000"])

        # Optional: minimum reads per cluster
        if "min_cluster_reads" in extra_args:
            cmd.extend(["-n", str(extra_args["min_cluster_reads"])])

        # Optional: ATL method
        atl_method = extra_args.get("atl_method", "p75")
        cmd.extend(["-m", atl_method])

        # minimap2 path override
        if "minimap2_path" in extra_args:
            cmd.extend(["--minimap2", extra_args["minimap2_path"]])

        returncode, stdout, stderr, elapsed = self._run_command(cmd, timeout=14400)

        result = ToolResult(
            tool_name=self.name,
            runtime_seconds=elapsed,
            raw_output_dir=output_dir,
        )

        if returncode != 0:
            result.error_message = stderr[:2000]
            return result

        result.success = True
        return self._parse_output(output_dir, result)

    def parse_results(self, output_dir: str) -> ToolResult:
        result = ToolResult(tool_name=self.name, raw_output_dir=output_dir)
        result.success = True
        return self._parse_output(output_dir, result)

    def _parse_output(self, output_dir: str, result: ToolResult) -> ToolResult:
        out = Path(output_dir)

        # --- Parse tlens_by_allele.tsv ---
        allele_tsv = out / "tlens_by_allele.tsv"
        if not allele_tsv.exists():
            result.error_message = f"Output not found: {allele_tsv}"
            result.success = False
            return result

        df = pd.read_csv(allele_tsv, sep="\t")

        # Identify the TL column (varies by -m flag: TL_p75, TL_mean, TL_median, TL_max)
        tl_col = None
        for col in df.columns:
            if col.startswith("TL_"):
                tl_col = col
                break
        if tl_col is None:
            result.error_message = "No TL column found in tlens_by_allele.tsv"
            result.success = False
            return result

        # Filter out interstitial alleles (allele_id contains 'i')
        df = df[~df["allele_id"].astype(str).str.contains("i", na=False)]
        # Filter out unmapped
        df = df[df["#chr"] != "chrU"]

        # --- Build per-read data by expanding comma-separated read_TLs ---
        read_rows = []
        for _, row in df.iterrows():
            arm = self._normalize_arm(str(row["#chr"]))
            allele_id = str(row["allele_id"])
            haplotype = self._infer_haplotype(str(row.get("ref_samp", "")))

            read_tls = str(row.get("read_TLs", ""))
            read_ids = str(row.get("supporting_reads", ""))

            tl_list = [self._safe_float(x) for x in read_tls.split(",") if x.strip()]
            id_list = [x.strip() for x in read_ids.split(",") if x.strip()]

            for j, tl in enumerate(tl_list):
                rid = id_list[j] if j < len(id_list) else f"{arm}_{allele_id}_read{j}"
                read_rows.append({
                    "read_id": rid,
                    "telomere_length": tl,
                    "chromosome_arm": arm,
                    "haplotype": haplotype,
                    "confidence": 1.0,  # Telogator2 doesn't output per-read confidence
                })

        if read_rows:
            result.read_level = pd.DataFrame(read_rows)

        # --- Build per-arm data ---
        arm_rows = []
        for _, row in df.iterrows():
            arm = self._normalize_arm(str(row["#chr"]))
            haplotype = self._infer_haplotype(str(row.get("ref_samp", "")))
            allele_id = str(row["allele_id"])

            read_tls = str(row.get("read_TLs", ""))
            tl_list = [self._safe_float(x) for x in read_tls.split(",") if x.strip()]

            if tl_list:
                arr = np.array(tl_list)
                arm_rows.append({
                    "chromosome_arm": arm,
                    "haplotype": haplotype,
                    "allele_id": allele_id,
                    "n_reads": len(arr),
                    "median_tl": float(np.median(arr)),
                    "mean_tl": float(np.mean(arr)),
                    "atl": self._safe_float(row[tl_col]),
                })

        if arm_rows:
            result.arm_level = pd.DataFrame(arm_rows)

        # --- Parse QC stats ---
        stats_tsv = out / "qc" / "stats.tsv"
        if stats_tsv.exists():
            stats = pd.read_csv(stats_tsv, sep="\t")
            if not stats.empty:
                s = stats.iloc[0]
                result.sample_summary = {
                    "median_tl": self._safe_float(s.get("tl_median")),
                    "mean_tl": self._safe_float(s.get("tl_mean")),
                    "n_reads": self._safe_int(s.get("num_telreads")),
                    "n_alleles": self._safe_int(s.get("num_alleles")),
                    "n_alleles_fail": self._safe_int(s.get("num_alleles_fail")),
                    "tl_short": self._safe_int(s.get("tl_short")),
                }
        elif result.arm_level is not None and not result.arm_level.empty:
            all_tl = result.read_level["telomere_length"] if result.read_level is not None else pd.Series()
            result.sample_summary = {
                "median_tl": float(all_tl.median()) if len(all_tl) > 0 else 0.0,
                "mean_tl": float(all_tl.mean()) if len(all_tl) > 0 else 0.0,
                "n_reads": len(all_tl),
                "n_arms": result.arm_level["chromosome_arm"].nunique(),
            }

        # --- Version ---
        rc, stdout, _, _ = self._run_command(
            ["telogator2", "--version"], timeout=10
        )
        if rc == 0 and stdout.strip():
            result.version = stdout.strip().split()[-1]

        return result

    @staticmethod
    def _infer_haplotype(ref_samp: str) -> str:
        """Infer haplotype from Telogator2 ref_samp field."""
        low = ref_samp.lower()
        if "pat" in low or "hap1" in low:
            return "pat"
        elif "mat" in low or "hap2" in low:
            return "mat"
        return "unknown"
