"""Wrapper for Topsicle (Nguyen & Choi, 2025).

Topsicle estimates telomere length from long reads using k-mer counting
and change-point detection via the ruptures library.

GitHub: https://github.com/jaeyoungchoilab/Topsicle
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .base import ToolWrapper, ToolResult

logger = logging.getLogger(__name__)


class TopsicleWrapper(ToolWrapper):
    """Wrapper for Topsicle."""

    def __init__(self):
        super().__init__("Topsicle")

    def check_installed(self) -> bool:
        return self._check_binary("topsicle")

    def run(
        self,
        input_path: str,
        output_dir: str,
        platform: str = "ont",
        reference: Optional[str] = None,
        threads: int = 4,
        extra_args: Optional[dict] = None,
    ) -> ToolResult:
        """Run Topsicle.

        Note: Topsicle requires FASTA/FASTQ input (not BAM).
              It does not perform chromosome-arm assignment.
              The --pattern must be the C-rich strand: CCCTAA for human.
        """
        extra_args = extra_args or {}
        pattern = extra_args.get("pattern", "CCCTAA")
        cutoff = extra_args.get("cutoff", 0.7)
        min_seq_length = extra_args.get("min_seq_length", 9000)
        window_size = extra_args.get("window_size", 100)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        cmd = [
            "topsicle",
            "--inputDir", input_path,
            "--outputDir", output_dir,
            "--pattern", pattern,
            "--threads", str(threads),
            "--cutoff", str(cutoff),
            "--minSeqLength", str(min_seq_length),
            "--windowSize", str(window_size),
        ]

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

        # --- Parse telolengths_all.csv ---
        telo_csv = out / "telolengths_all.csv"
        if not telo_csv.exists():
            result.error_message = f"Output not found: {telo_csv}"
            result.success = False
            return result

        df = pd.read_csv(telo_csv)

        # Expected columns: file_number, phrase, trc, readID, telo_length
        if "telo_length" not in df.columns or "readID" not in df.columns:
            result.error_message = f"Unexpected columns in telolengths_all.csv: {list(df.columns)}"
            result.success = False
            return result

        # Filter valid telomere lengths
        df = df[df["telo_length"] > 0].copy()

        # --- Build per-read data ---
        # Topsicle doesn't do arm assignment, so chromosome_arm is empty
        result.read_level = pd.DataFrame({
            "read_id": df["readID"],
            "telomere_length": df["telo_length"].astype(float),
            "chromosome_arm": "",
            "haplotype": "",
            "confidence": df["trc"].astype(float),  # TRC as confidence proxy
        })

        # --- Sample summary ---
        tl = df["telo_length"].astype(float)
        result.sample_summary = {
            "median_tl": float(tl.median()),
            "mean_tl": float(tl.mean()),
            "n_reads": len(tl),
            "n_arms": 0,  # Topsicle doesn't assign arms
            "trc_median": float(df["trc"].median()),
            "trc_mean": float(df["trc"].mean()),
        }

        # No arm-level data (Topsicle doesn't perform arm assignment)
        result.arm_level = None

        # --- Version ---
        rc, stdout, _, _ = self._run_command(["topsicle", "--version"], timeout=10)
        if rc == 0 and stdout.strip():
            result.version = stdout.strip()

        # --- Parse log for additional metrics ---
        log_file = out / "topsicle_run.log"
        if log_file.exists():
            self._parse_log(log_file, result)

        return result

    def _parse_log(self, log_path: Path, result: ToolResult) -> None:
        """Extract additional metrics from Topsicle's run log."""
        try:
            text = log_path.read_text()
            for line in text.split("\n"):
                if "Median telomere length" in line:
                    parts = line.split(":")
                    if len(parts) >= 2:
                        result.sample_summary["log_median_tl"] = self._safe_float(
                            parts[-1].strip().replace("bp", "").strip()
                        )
                if "asymptotic TRC" in line.lower():
                    parts = line.split(":")
                    if len(parts) >= 2:
                        result.sample_summary["asymptotic_trc"] = self._safe_float(
                            parts[-1].strip()
                        )
        except Exception as e:
            logger.debug("Failed to parse Topsicle log: %s", e)
