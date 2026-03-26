"""Sample sheet parsing and validation."""

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ..constants import Mode, Platform

logger = logging.getLogger(__name__)


@dataclass
class SampleConfig:
    """Configuration for a single sample."""

    sample_id: str
    input_path: str
    barcode: Optional[str] = None
    reference: Optional[str] = None
    mode: Mode = Mode.TELOSEQ
    platform: Platform = Platform.ONT


REQUIRED_COLUMNS = {"sample_id", "input_path"}
OPTIONAL_COLUMNS = {"barcode", "reference", "mode", "platform"}


def parse_sample_sheet(path: str | Path) -> list[SampleConfig]:
    """Parse a CSV sample sheet into SampleConfig objects."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Sample sheet not found: {path}")

    samples = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        columns = set(reader.fieldnames or [])

        missing = REQUIRED_COLUMNS - columns
        if missing:
            raise ValueError(
                f"Sample sheet missing required columns: {missing}"
            )

        for i, row in enumerate(reader, start=2):
            sample_id = row["sample_id"].strip()
            if not sample_id:
                raise ValueError(f"Empty sample_id at row {i}")

            mode_str = row.get("mode", "teloprobe").strip().lower()
            platform_str = row.get("platform", "ont").strip().lower()

            samples.append(SampleConfig(
                sample_id=sample_id,
                input_path=row["input_path"].strip(),
                barcode=row.get("barcode", "").strip() or None,
                reference=row.get("reference", "").strip() or None,
                mode=Mode(mode_str) if mode_str else Mode.TELOSEQ,
                platform=Platform(platform_str) if platform_str else Platform.ONT,
            ))

    logger.info("Parsed %d samples from %s", len(samples), path)
    return samples


def validate_sample_sheet(samples: list[SampleConfig]) -> list[str]:
    """Validate sample sheet entries. Returns list of error messages."""
    errors = []
    seen_ids = set()

    for s in samples:
        if s.sample_id in seen_ids:
            errors.append(f"Duplicate sample_id: {s.sample_id}")
        seen_ids.add(s.sample_id)

        if not Path(s.input_path).exists():
            errors.append(
                f"Input path not found for {s.sample_id}: {s.input_path}"
            )

        if s.reference and not Path(s.reference).exists():
            errors.append(
                f"Reference not found for {s.sample_id}: {s.reference}"
            )

    return errors
