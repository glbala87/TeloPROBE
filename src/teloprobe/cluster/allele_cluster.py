"""Allele-specific telomere clustering using subtelomeric flanks and TVR profiles.

Uses GMM or HDBSCAN to separate reads from the same chromosome arm
into allele groups based on:
- Telomere length distribution (bimodality)
- TVR composition profile
- Subtelomeric flank similarity
"""

import logging
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def cluster_alleles(
    tl_values: np.ndarray,
    tvr_features: Optional[np.ndarray] = None,
    method: str = "gmm",
    n_components: int = 2,
    min_reads: int = 10,
) -> np.ndarray:
    """Cluster reads into allele groups.

    Args:
        tl_values: Telomere length values for reads at a single arm.
        tvr_features: Optional TVR profile features (n_reads x n_features).
        method: Clustering method - "gmm" or "hdbscan".
        n_components: Number of expected alleles (default 2 for diploid).
        min_reads: Minimum reads required for clustering.

    Returns:
        Array of cluster labels (0, 1, ...) for each read.
        Returns all-zeros if clustering is not justified.
    """
    n = len(tl_values)
    if n < min_reads:
        return np.zeros(n, dtype=int)

    # Build feature matrix
    features = tl_values.reshape(-1, 1)
    if tvr_features is not None and tvr_features.shape[0] == n:
        # Normalize TVR features to same scale as TL
        from sklearn.preprocessing import StandardScaler
        tvr_scaled = StandardScaler().fit_transform(tvr_features)
        tl_scaled = StandardScaler().fit_transform(features)
        features = np.hstack([tl_scaled, tvr_scaled * 0.5])
    else:
        from sklearn.preprocessing import StandardScaler
        features = StandardScaler().fit_transform(features)

    if method == "gmm":
        return _cluster_gmm(features, n_components)
    elif method == "hdbscan":
        return _cluster_hdbscan(features)
    else:
        logger.warning("Unknown clustering method '%s', using gmm", method)
        return _cluster_gmm(features, n_components)


def _cluster_gmm(features: np.ndarray, n_components: int = 2) -> np.ndarray:
    """Gaussian Mixture Model clustering."""
    from sklearn.mixture import GaussianMixture

    # Compare 1-component vs 2-component via BIC
    gmm1 = GaussianMixture(n_components=1, random_state=42).fit(features)
    gmm2 = GaussianMixture(n_components=n_components, random_state=42).fit(features)

    # Only use 2 components if BIC improves substantially
    if gmm2.bic(features) < gmm1.bic(features) - 10:
        labels = gmm2.predict(features)
    else:
        labels = np.zeros(len(features), dtype=int)

    return labels


def _cluster_hdbscan(features: np.ndarray) -> np.ndarray:
    """HDBSCAN clustering (density-based, no fixed k)."""
    try:
        from sklearn.cluster import HDBSCAN
        clusterer = HDBSCAN(min_cluster_size=5, min_samples=3)
        labels = clusterer.fit_predict(features)
        # Convert noise (-1) to 0
        labels[labels == -1] = 0
        return labels
    except ImportError:
        logger.warning("HDBSCAN not available, falling back to GMM")
        return _cluster_gmm(features, 2)


def is_bimodal(
    tl_values: np.ndarray,
    min_reads: int = 20,
) -> bool:
    """Test if telomere length distribution is bimodal.

    Uses BIC comparison between 1-component and 2-component GMMs.
    """
    if len(tl_values) < min_reads:
        return False

    from sklearn.mixture import GaussianMixture

    X = tl_values.reshape(-1, 1)
    gmm1 = GaussianMixture(n_components=1, random_state=42).fit(X)
    gmm2 = GaussianMixture(n_components=2, random_state=42).fit(X)

    # BIC improvement threshold
    return gmm2.bic(X) < gmm1.bic(X) - 20


def evaluate_clustering(
    tl_values: np.ndarray,
    labels: np.ndarray,
) -> float:
    """Evaluate clustering quality using silhouette score.

    Returns score between -1 and 1 (higher = better separation).
    Returns 0.0 if only one cluster.
    """
    n_clusters = len(set(labels))
    if n_clusters <= 1 or n_clusters >= len(tl_values):
        return 0.0

    from sklearn.metrics import silhouette_score
    return silhouette_score(tl_values.reshape(-1, 1), labels)


def cluster_arms(
    read_level: pd.DataFrame,
    method: str = "gmm",
    min_reads: int = 10,
) -> pd.DataFrame:
    """Cluster alleles for all chromosome arms in the read-level data.

    Returns DataFrame with allele_cluster column added.
    """
    result = read_level.copy()
    result["allele_cluster"] = 0

    arms = result[result["chromosome_arm"].notna()]["chromosome_arm"].unique()

    for arm in arms:
        mask = result["chromosome_arm"] == arm
        arm_data = result[mask]

        if len(arm_data) < min_reads:
            continue

        tl = arm_data["telomere_length_bp"].values.astype(float)
        valid = tl > 0
        if valid.sum() < min_reads:
            continue

        labels = cluster_alleles(tl[valid], method=method, min_reads=min_reads)

        # Map labels back
        valid_idx = arm_data.index[valid]
        for idx, label in zip(valid_idx, labels):
            result.loc[idx, "allele_cluster"] = int(label)

    return result
