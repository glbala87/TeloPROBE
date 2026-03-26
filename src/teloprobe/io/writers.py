"""TSV and JSON output writers."""

import json
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def write_read_level(df: pd.DataFrame, path: str | Path) -> None:
    """Write read-level results to TSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    logger.info("Wrote read-level results (%d rows) to %s", len(df), path)


def write_arm_level(df: pd.DataFrame, path: str | Path) -> None:
    """Write arm-level results to TSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    logger.info("Wrote arm-level results (%d rows) to %s", len(df), path)


def write_sample_level(df: pd.DataFrame, path: str | Path) -> None:
    """Write sample-level results to TSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    logger.info("Wrote sample-level results to %s", path)


def write_qc_json(qc_checks: list, path: str | Path) -> None:
    """Write QC check results to JSON."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    data = []
    for check in qc_checks:
        data.append({
            "name": check.name,
            "passed": check.passed,
            "value": float(check.value),
            "threshold": float(check.threshold),
            "message": check.message,
        })

    path.write_text(json.dumps(data, indent=2))
    logger.info("Wrote QC checks to %s", path)


def write_config_json(config, path: str | Path) -> None:
    """Write pipeline configuration to JSON for reproducibility."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(config.to_dict(), indent=2, default=str))
    logger.info("Wrote config to %s", path)
