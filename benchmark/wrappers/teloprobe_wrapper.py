"""Wrapper for TeloPROBE (this tool)."""

import logging
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .base import ToolWrapper, ToolResult

logger = logging.getLogger(__name__)


class TeloPROBEWrapper(ToolWrapper):
    """Wrapper for running TeloPROBE in benchmark mode."""

    def __init__(self):
        super().__init__("TeloPROBE")

    def check_installed(self) -> bool:
        return self._check_binary("teloprobe")

    def run(
        self,
        input_path: str,
        output_dir: str,
        platform: str = "ont",
        reference: Optional[str] = None,
        threads: int = 4,
        extra_args: Optional[dict] = None,
    ) -> ToolResult:
        """Run TeloPROBE via CLI."""
        extra_args = extra_args or {}
        mode = extra_args.get("mode", "teloprobe")
        sample_id = extra_args.get("sample_id", "benchmark")

        cmd = [
            "teloprobe", "run",
            "-i", input_path,
            "-o", output_dir,
            "--platform", platform,
            "--mode", mode,
            "--sample-id", sample_id,
            "-t", str(threads),
        ]

        if reference:
            cmd.extend(["-r", reference])
        else:
            cmd.append("--skip-mapping")

        for key, val in extra_args.items():
            if key not in ("mode", "sample_id"):
                cmd.extend([f"--{key}", str(val)])

        returncode, stdout, stderr, elapsed = self._run_command(cmd)

        result = ToolResult(
            tool_name=self.name,
            runtime_seconds=elapsed,
            raw_output_dir=output_dir,
        )

        if returncode != 0:
            result.error_message = stderr[:2000]
            return result

        result.success = True
        return self._parse_output(output_dir, sample_id, result)

    def parse_results(self, output_dir: str) -> ToolResult:
        result = ToolResult(tool_name=self.name, raw_output_dir=output_dir)
        # Find sample subdirectory
        out = Path(output_dir)
        sample_dirs = [d for d in out.iterdir() if d.is_dir() and (d / "read_level.tsv").exists()]
        if not sample_dirs:
            result.error_message = "No sample output found"
            return result
        sample_id = sample_dirs[0].name
        result.success = True
        return self._parse_output(output_dir, sample_id, result)

    def _parse_output(self, output_dir: str, sample_id: str, result: ToolResult) -> ToolResult:
        base = Path(output_dir) / sample_id

        # Read-level
        rl_path = base / "read_level.tsv"
        if rl_path.exists():
            df = pd.read_csv(rl_path, sep="\t")
            good = df[df["qc_status"] == "Good"].copy()
            result.read_level = pd.DataFrame({
                "read_id": good["read_id"],
                "telomere_length": good["telomere_length_bp"],
                "chromosome_arm": good.get("chromosome_arm", ""),
                "haplotype": good.get("haplotype", ""),
                "confidence": good.get("confidence_score", 0.0),
            })

        # Arm-level
        al_path = base / "arm_level.tsv"
        if al_path.exists() and al_path.stat().st_size > 10:
            try:
                df = pd.read_csv(al_path, sep="\t")
            except Exception:
                df = pd.DataFrame()
            if not df.empty:
                result.arm_level = df.rename(columns={
                    "median_tl": "median_tl",
                    "mean_tl": "mean_tl",
                })

        # Sample-level
        sl_path = base / "sample_level.tsv"
        if sl_path.exists():
            df = pd.read_csv(sl_path, sep="\t")
            if not df.empty:
                row = df.iloc[0]
                result.sample_summary = {
                    "median_tl": self._safe_float(row.get("median_tl")),
                    "mean_tl": self._safe_float(row.get("trimmed_mean_tl")),
                    "n_reads": self._safe_int(row.get("informative_reads")),
                    "total_reads": self._safe_int(row.get("total_reads")),
                    "n_arms": self._safe_int(row.get("arms_assigned")),
                    "ci_lower": self._safe_float(row.get("global_ci_lower")),
                    "ci_upper": self._safe_float(row.get("global_ci_upper")),
                    "cv": self._safe_float(row.get("cv")),
                }

        # Version
        try:
            from teloprobe import __version__
            result.version = __version__
        except ImportError:
            result.version = "unknown"

        return result
