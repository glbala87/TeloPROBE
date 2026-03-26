"""Platform detection from BAM headers and read characteristics."""

import logging
from pathlib import Path

import pysam

from ..constants import Platform, PLATFORM_PRESETS

logger = logging.getLogger(__name__)

# Keywords in BAM headers that identify platforms
ONT_KEYWORDS = {"dorado", "guppy", "minknow", "nanopore", "ont"}
PACBIO_KEYWORDS = {"ccs", "activ", "sequel", "revio", "pacbio", "smrtlink", "ics"}


def detect_platform(path: str | Path) -> Platform:
    """Detect sequencing platform from BAM/CRAM header.

    Inspects @RG (read group) and @PG (program) header lines for
    platform-specific keywords.
    """
    path = Path(path)
    if not path.exists():
        return Platform.UNKNOWN

    try:
        with pysam.AlignmentFile(str(path), "rb", check_sq=False) as af:
            header = af.header.to_dict()
    except Exception:
        return Platform.UNKNOWN

    # Check read groups for platform (PL) tag
    for rg in header.get("RG", []):
        pl = rg.get("PL", "").lower()
        if any(kw in pl for kw in ONT_KEYWORDS):
            logger.info("Detected ONT platform from RG:PL=%s", rg.get("PL"))
            return Platform.ONT
        if any(kw in pl for kw in PACBIO_KEYWORDS):
            logger.info("Detected PacBio platform from RG:PL=%s", rg.get("PL"))
            return Platform.PACBIO_HIFI

    # Check program lines
    for pg in header.get("PG", []):
        pn = pg.get("PN", "").lower()
        cl = pg.get("CL", "").lower()
        combined = f"{pn} {cl}"
        if any(kw in combined for kw in ONT_KEYWORDS):
            logger.info("Detected ONT platform from PG header")
            return Platform.ONT
        if any(kw in combined for kw in PACBIO_KEYWORDS):
            logger.info("Detected PacBio platform from PG header")
            return Platform.PACBIO_HIFI

    logger.warning("Could not detect platform from %s, defaulting to UNKNOWN", path)
    return Platform.UNKNOWN


def get_platform_preset(platform: Platform) -> dict:
    """Return platform-specific parameter preset."""
    return dict(PLATFORM_PRESETS.get(platform, PLATFORM_PRESETS[Platform.ONT]))
