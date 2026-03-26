"""Mixture model fitting for telomere length distributions.

Fits Gaussian mixture models to telomere length distributions
with automatic component selection via BIC.
"""

import logging
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class MixtureResult:
    """Result of mixture model fitting."""

    n_components: int
    means: np.ndarray
    stds: np.ndarray
    weights: np.ndarray
    bic: float
    converged: bool = True

    @property
    def dominant_mean(self) -> float:
        """Mean of the component with the highest weight."""
        idx = np.argmax(self.weights)
        return float(self.means[idx])


def fit_mixture(
    values: np.ndarray,
    max_components: int = 3,
    min_reads: int = 20,
) -> Optional[MixtureResult]:
    """Fit a Gaussian mixture model with BIC-based component selection.

    Args:
        values: Telomere length values.
        max_components: Maximum number of components to try.
        min_reads: Minimum reads required.

    Returns:
        MixtureResult with the best-fitting model, or None if insufficient data.
    """
    if len(values) < min_reads:
        return None

    from sklearn.mixture import GaussianMixture

    X = values.reshape(-1, 1)
    best_result = None
    best_bic = np.inf

    for k in range(1, max_components + 1):
        if k > len(values) // 5:
            break

        gmm = GaussianMixture(
            n_components=k,
            covariance_type="full",
            random_state=42,
            n_init=3,
        ).fit(X)

        bic = gmm.bic(X)

        if bic < best_bic:
            best_bic = bic
            best_result = MixtureResult(
                n_components=k,
                means=gmm.means_.flatten(),
                stds=np.sqrt(gmm.covariances_.flatten()),
                weights=gmm.weights_,
                bic=bic,
                converged=gmm.converged_,
            )

    return best_result


def assign_components(
    values: np.ndarray,
    result: MixtureResult,
) -> np.ndarray:
    """Assign values to mixture components using posterior probabilities."""
    from sklearn.mixture import GaussianMixture

    X = values.reshape(-1, 1)
    gmm = GaussianMixture(
        n_components=result.n_components,
        random_state=42,
    )
    gmm.means_ = result.means.reshape(-1, 1)
    gmm.covariances_ = (result.stds ** 2).reshape(-1, 1, 1)
    gmm.weights_ = result.weights
    gmm.precisions_cholesky_ = np.sqrt(
        1.0 / gmm.covariances_
    ).reshape(-1, 1, 1)

    return gmm.predict(X)
