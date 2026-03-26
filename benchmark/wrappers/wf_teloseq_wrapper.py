"""Wrapper for wf-teloseq (EPI2ME Labs / ONT, Nextflow pipeline).

wf-teloseq is a Nextflow-based pipeline for TeloPROBE-enriched ONT data.
It performs telomere boundary detection and optional alignment to a
subtelomeric reference.

GitHub: https://github.com/epi2me-labs/wf-teloseq
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .base import ToolWrapper, ToolResult

logger = logging.getLogger(__name__)


class WfTeloseqWrapper(ToolWrapper):
    """Wrapper for wf-teloseq (Nextflow pipeline)."""

    def __init__(self):
        super().__init__("wf-teloseq")

    def check_installed(self) -> bool:
        """Check if Nextflow is available."""
        return self._check_binary("nextflow")

    def run(
        self,
        input_path: str,
        output_dir: str,
        platform: str = "ont",
        reference: Optional[str] = None,
        threads: int = 4,
        extra_args: Optional[dict] = None,
    ) -> ToolResult:
        """Run wf-teloseq via Nextflow.

        Args:
            input_path: BAM or FASTQ input file or directory.
            reference: Subtelomeric reference FASTA (optional).
            extra_args: May include 'workflow_dir' pointing to the wf-teloseq
                        installation, 'sample' name, 'profile' (docker/local).
        """
        extra_args = extra_args or {}
        workflow_dir = extra_args.get(
            "workflow_dir",
            str(Path(__file__).parents[2].parent)  # Default: parent of teloprobe
        )
        sample_name = extra_args.get("sample", "benchmark")
        profile = extra_args.get("profile", "standard")

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Determine if input is BAM or FASTQ
        input_path_p = Path(input_path)
        if input_path_p.suffix in (".bam", ".cram"):
            input_flag = "--bam"
        else:
            input_flag = "--fastq"

        cmd = [
            "nextflow", "run", workflow_dir,
            input_flag, input_path,
            "--out_dir", output_dir,
            "--sample", sample_name,
            "-profile", profile,
        ]

        if reference:
            cmd.extend(["--reference", reference])
        else:
            cmd.append("--skip_mapping")

        # Thread control
        cmd.extend(["--alignment_threads", str(threads)])

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
        return self._parse_output(output_dir, sample_name, result)

    def parse_results(self, output_dir: str) -> ToolResult:
        result = ToolResult(tool_name=self.name, raw_output_dir=output_dir)
        result.success = True
        # Find sample name from directory structure
        out = Path(output_dir)
        sample_dirs = [
            d for d in out.iterdir()
            if d.is_dir() and (d / "stats").is_dir()
        ]
        if sample_dirs:
            return self._parse_output(output_dir, sample_dirs[0].name, result)
        result.error_message = "No sample directory found"
        result.success = False
        return result

    def _parse_output(self, output_dir: str, sample: str, result: ToolResult) -> ToolResult:
        out = Path(output_dir) / sample
        stats_dir = out / "stats"

        # --- Parse per-read from BAM/FASTQ tags ---
        # wf-teloseq stores TL in tl:i tags in the output FASTQ
        # We parse the unaligned metrics instead for summary
        unaligned = stats_dir / f"{sample}_telomere_unaligned_metrics.tsv"
        if unaligned.exists():
            df = pd.read_csv(unaligned, sep="\t")
            if not df.empty:
                row = df.iloc[0]
                result.sample_summary = {
                    "median_tl": self._safe_float(row.get("Median length")),
                    "mean_tl": 0.0,  # Not directly reported
                    "n_reads": self._safe_int(row.get("Read count")),
                    "cv": self._safe_float(row.get("CV")),
                    "q1": self._safe_float(row.get("Q1")),
                    "q3": self._safe_float(row.get("Q3")),
                }

        # --- Parse per-contig (arm-level) ---
        contig_file = stats_dir / f"{sample}_contig_telomere_aligned_metrics.tsv"
        if contig_file.exists():
            df = pd.read_csv(contig_file, sep="\t")
            if not df.empty:
                arm_rows = []
                for _, row in df.iterrows():
                    contig = str(row.iloc[0])  # First column is contig name
                    arm = self._parse_contig_name(contig)
                    haplotype = self._parse_haplotype(contig)

                    arm_rows.append({
                        "chromosome_arm": arm,
                        "haplotype": haplotype,
                        "allele_id": "",
                        "n_reads": self._safe_int(row.get("Read count")),
                        "median_tl": self._safe_float(row.get("Median length")),
                        "mean_tl": 0.0,
                    })

                if arm_rows:
                    result.arm_level = pd.DataFrame(arm_rows)
                    result.sample_summary["n_arms"] = len(arm_rows)

        # --- Parse per-read from output FASTQ tags ---
        fastq_path = out / "unaligned_data" / f"{sample}_filtered_telomeric.fastq"
        if fastq_path.exists():
            result.read_level = self._parse_tagged_fastq(fastq_path)

        # --- Parse QC modes ---
        qc_file = stats_dir / f"{sample}_qc_modes_metrics.tsv"
        if qc_file.exists():
            df = pd.read_csv(qc_file, sep="\t")
            if not df.empty:
                result.sample_summary["qc_breakdown"] = dict(
                    zip(df.iloc[:, 0], df.iloc[:, 1])
                )

        result.version = "1.0.4"
        return result

    def _parse_tagged_fastq(self, fastq_path: Path) -> Optional[pd.DataFrame]:
        """Parse SAM-tagged FASTQ to extract per-read TL and QC status."""
        try:
            import pysam
        except ImportError:
            logger.warning("pysam not available for FASTQ parsing")
            return None

        rows = []
        try:
            with pysam.FastxFile(str(fastq_path)) as fq:
                for record in fq:
                    comment = record.comment or ""
                    tl = -1
                    qc = "Unknown"

                    # Parse SAM tags from comment: tl:i:1234 qc:Z:Good
                    for tag in comment.split("\t"):
                        if tag.startswith("tl:i:"):
                            tl = self._safe_int(tag[5:], -1)
                        elif tag.startswith("qc:Z:"):
                            qc = tag[5:]

                    if qc == "Good" and tl > 0:
                        rows.append({
                            "read_id": record.name,
                            "telomere_length": float(tl),
                            "chromosome_arm": "",
                            "haplotype": "",
                            "confidence": 1.0,
                        })
        except Exception as e:
            logger.warning("Failed to parse tagged FASTQ: %s", e)
            return None

        if rows:
            return pd.DataFrame(rows)
        return None

    @staticmethod
    def _parse_contig_name(contig: str) -> str:
        """Parse wf-teloseq contig name to chromosome arm.
        E.g.: chr1_mat_p -> chr1p, chr1_PATERNAL_P -> chr1p
        """
        import re
        m = re.match(
            r"(chr\d+|chrX|chrY)[-_](?:pat(?:ernal)?|mat(?:ernal)?|hap[12])[-_]([pq])",
            contig, re.IGNORECASE
        )
        if m:
            return f"{m.group(1).lower()}{m.group(2).lower()}"

        m = re.match(r"(chr\d+|chrX|chrY)([pq])", contig, re.IGNORECASE)
        if m:
            return f"{m.group(1).lower()}{m.group(2).lower()}"

        return contig.lower()

    @staticmethod
    def _parse_haplotype(contig: str) -> str:
        low = contig.lower()
        if "pat" in low or "hap1" in low:
            return "pat"
        elif "mat" in low or "hap2" in low:
            return "mat"
        return "unknown"
