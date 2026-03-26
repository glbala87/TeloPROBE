"""Configuration dataclasses and defaults."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml

from .constants import Mode, Platform, PLATFORM_PRESETS


@dataclass
class Config:
    """Main configuration for teloprobe pipeline."""

    # Input/output
    input_path: Optional[str] = None
    output_dir: str = "output"
    reference: Optional[str] = None
    sample_sheet: Optional[str] = None
    sample_id: str = "sample"

    # Mode and platform
    mode: Mode = Mode.TELOSEQ
    platform: Platform = Platform.ONT

    # Processing
    threads: int = 4
    skip_mapping: bool = False
    random_seed: int = 42

    # Read filtering
    min_quality: int = 9
    min_length: int = 100
    min_read_length_boundary: int = 160

    # Motif scanning
    min_repeats: int = 100
    filter_width: int = 10
    motif_density_threshold: float = 0.6
    motif_window_size: int = 120

    # Boundary detection
    start_window_frac: float = 0.3
    start_repeats_frac: float = 0.8
    min_qual_non_telo: int = 9
    post_boundary_ccc_threshold: float = 0.25
    max_errors: int = 5
    error_distance: int = 500

    # Change-point detection
    changepoint_method: str = "pelt"
    changepoint_penalty: Optional[float] = None
    changepoint_min_size: int = 60

    # HMM
    hmm_n_states: int = 3
    hmm_window_size: int = 6

    # Alignment
    minimap2_preset: str = "lr:hq"
    identity_threshold: float = 0.8
    mapq_threshold: int = 20

    # Clustering
    allele_min_reads: int = 10
    allele_method: str = "gmm"

    # Confidence / bootstrap
    bootstrap_n: int = 1000
    bootstrap_alpha: float = 0.05

    # QC
    min_informative_reads: int = 20
    min_arm_reads: int = 10

    # Report
    offline: bool = False

    # Valid options for string-valued parameters
    VALID_CHANGEPOINT_METHODS = {"pelt", "binseg", "window", "dynp"}
    VALID_ALLELE_METHODS = {"gmm", "kmeans", "manual"}

    def validate(self):
        """Validate configuration parameters.

        Raises:
            ValueError: If any parameter is out of its valid range.
        """
        if self.bootstrap_n <= 0:
            raise ValueError(
                f"bootstrap_n must be > 0, got {self.bootstrap_n}"
            )
        if self.min_length <= 0:
            raise ValueError(
                f"min_length must be > 0, got {self.min_length}"
            )
        if self.min_repeats <= 0:
            raise ValueError(
                f"min_repeats must be > 0, got {self.min_repeats}"
            )
        if self.threads <= 0:
            raise ValueError(
                f"threads must be > 0, got {self.threads}"
            )
        if not (0 < self.bootstrap_alpha < 1):
            raise ValueError(
                f"bootstrap_alpha must be between 0 and 1 (exclusive), "
                f"got {self.bootstrap_alpha}"
            )
        if self.changepoint_method not in self.VALID_CHANGEPOINT_METHODS:
            raise ValueError(
                f"changepoint_method must be one of "
                f"{self.VALID_CHANGEPOINT_METHODS}, "
                f"got '{self.changepoint_method}'"
            )
        if self.allele_method not in self.VALID_ALLELE_METHODS:
            raise ValueError(
                f"allele_method must be one of "
                f"{self.VALID_ALLELE_METHODS}, "
                f"got '{self.allele_method}'"
            )

    def apply_platform_preset(self):
        """Override defaults with platform-specific presets."""
        preset = PLATFORM_PRESETS.get(self.platform, {})
        for key, value in preset.items():
            if hasattr(self, key):
                setattr(self, key, value)
        self.validate()

    @classmethod
    def from_yaml(cls, path: str | Path) -> "Config":
        """Load configuration from YAML file."""
        with open(path) as f:
            data = yaml.safe_load(f) or {}

        # Convert string enums
        if "mode" in data:
            data["mode"] = Mode(data["mode"])
        if "platform" in data:
            data["platform"] = Platform(data["platform"])

        valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in data.items() if k in valid_fields}
        config = cls(**filtered)
        config.apply_platform_preset()
        config.validate()
        return config

    def to_dict(self) -> dict:
        """Serialize to dict for Snakemake config."""
        result = {}
        for k, v in self.__dict__.items():
            if isinstance(v, (Mode, Platform)):
                result[k] = v.value
            else:
                result[k] = v
        return result
