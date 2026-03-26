"""Base class for tool wrappers."""

import logging
import shutil
import subprocess
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class ToolResult:
    """Standardized result from any telomere tool."""

    tool_name: str
    version: str = ""
    runtime_seconds: float = 0.0
    peak_memory_mb: float = 0.0
    success: bool = False
    error_message: str = ""

    # Standardized per-read data
    read_level: Optional[pd.DataFrame] = None
    # Columns: read_id, telomere_length, chromosome_arm, haplotype, confidence

    # Standardized per-arm data
    arm_level: Optional[pd.DataFrame] = None
    # Columns: chromosome_arm, haplotype, allele_id, n_reads, median_tl, mean_tl

    # Standardized sample summary
    sample_summary: dict = field(default_factory=dict)
    # Keys: median_tl, mean_tl, n_reads, n_arms, ...

    raw_output_dir: Optional[str] = None


class ToolWrapper(ABC):
    """Abstract base class for tool wrappers.

    Each wrapper must implement:
      - check_installed(): verify the tool is available
      - run(): execute the tool on input data
      - parse_results(): parse outputs into standardized ToolResult
    """

    def __init__(self, name: str):
        self.name = name

    @abstractmethod
    def check_installed(self) -> bool:
        """Check if the tool is installed and available."""
        ...

    @abstractmethod
    def run(
        self,
        input_path: str,
        output_dir: str,
        platform: str = "ont",
        reference: Optional[str] = None,
        threads: int = 4,
        extra_args: Optional[dict] = None,
    ) -> ToolResult:
        """Run the tool and return standardized results."""
        ...

    @abstractmethod
    def parse_results(self, output_dir: str) -> ToolResult:
        """Parse existing outputs into standardized ToolResult."""
        ...

    def _run_command(
        self,
        cmd: list[str],
        timeout: int = 7200,
        capture_output: bool = True,
    ) -> tuple[int, str, str, float]:
        """Run a shell command, returning (returncode, stdout, stderr, runtime_seconds)."""
        logger.info("[%s] Running: %s", self.name, " ".join(cmd))
        start = time.time()

        try:
            result = subprocess.run(
                cmd,
                capture_output=capture_output,
                text=True,
                timeout=timeout,
            )
            elapsed = time.time() - start
            return result.returncode, result.stdout, result.stderr, elapsed

        except subprocess.TimeoutExpired:
            elapsed = time.time() - start
            return -1, "", f"Timeout after {timeout}s", elapsed

        except FileNotFoundError as e:
            elapsed = time.time() - start
            return -1, "", str(e), elapsed

    def _check_binary(self, binary: str) -> bool:
        """Check if a binary is on PATH."""
        return shutil.which(binary) is not None

    @staticmethod
    def _safe_float(val, default: float = 0.0) -> float:
        try:
            return float(val)
        except (ValueError, TypeError):
            return default

    @staticmethod
    def _safe_int(val, default: int = 0) -> int:
        try:
            return int(float(val))
        except (ValueError, TypeError):
            return default

    @staticmethod
    def _normalize_arm(arm_str: str) -> str:
        """Normalize chromosome arm names to consistent format: chr1p, chr1q, etc."""
        if not arm_str or arm_str in ("chrU", ""):
            return ""
        arm = arm_str.strip().lower()
        # Remove _ separators: chr1_p -> chr1p
        arm = arm.replace("_", "")
        return arm
