"""HMM-based window labeler for telomere/subtelomere classification.

Uses a 3-state Gaussian HMM:
  State 0: Telomeric (high motif density ~0.85-1.0)
  State 1: Transition zone (medium density ~0.3-0.7)
  State 2: Subtelomeric/genomic (low density ~0.0-0.15)

The HMM provides a generative model that handles sequencing error
and variant repeat structure better than hard thresholds.
"""

import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

try:
    from hmmlearn.hmm import GaussianHMM
    HAS_HMMLEARN = True
except ImportError:
    HAS_HMMLEARN = False
    logger.warning("hmmlearn not installed; HMM labeling disabled")


def build_hmm(
    n_states: int = 3,
    telo_mean: float = 0.85,
    subtelo_mean: float = 0.05,
    random_seed: int = 42,
) -> Optional["GaussianHMM"]:
    """Create a 3-state HMM with biologically informed fixed priors.

    This HMM uses fixed, biologically-informed parameters rather than
    data-driven training.  We intentionally set ``init_params=""`` and
    ``params=""`` so that hmmlearn will not attempt to re-estimate any
    parameters during decoding.  This is the correct approach because:

    1. Individual telomeric reads are too short for reliable parameter
       estimation -- a single read typically covers only one transition
       from telomeric to subtelomeric sequence.
    2. The emission and transition distributions are well-characterised
       from prior biological knowledge of telomere repeat structure.
    3. Keeping parameters fixed ensures reproducibility across samples
       and avoids overfitting to per-read noise.

    If data-driven parameter estimation is desired (e.g. on a large batch
    of density vectors), use :func:`fit_hmm` instead.

    The transition matrix encodes that:
    - Telomeric state is sticky (high self-transition)
    - Transition state is brief
    - Subtelomeric state is sticky
    - Transitions go telo -> transition -> subtelo (left to right bias)

    Args:
        n_states: Number of HMM states (default 3).
        telo_mean: Expected motif density in telomeric regions.
        subtelo_mean: Expected motif density in subtelomeric regions.
        random_seed: For reproducibility.

    Returns:
        Initialized GaussianHMM, or None if hmmlearn not available.
    """
    if not HAS_HMMLEARN:
        return None

    model = GaussianHMM(
        n_components=n_states,
        covariance_type="diag",
        n_iter=100,
        random_state=random_seed,
        init_params="",  # We set all params manually below
        params="",       # Do NOT train -- use fixed biologically-informed priors
    )

    # Start probabilities: reads start in telomeric state
    model.startprob_ = np.array([0.9, 0.05, 0.05])

    # Transition matrix: left-to-right bias
    # telo -> telo: 0.95, telo -> transition: 0.04, telo -> subtelo: 0.01
    # transition -> telo: 0.05, transition -> transition: 0.3, transition -> subtelo: 0.65
    # subtelo -> telo: 0.01, subtelo -> transition: 0.02, subtelo -> subtelo: 0.97
    model.transmat_ = np.array([
        [0.95, 0.04, 0.01],
        [0.05, 0.30, 0.65],
        [0.01, 0.02, 0.97],
    ])

    # Emission means (motif density per state)
    transition_mean = (telo_mean + subtelo_mean) / 2
    model.means_ = np.array([
        [telo_mean],
        [transition_mean],
        [subtelo_mean],
    ])

    # Emission variances (shape: n_components x n_features for 'diag')
    model.covars_ = np.array([
        [0.02],   # telomeric: tight distribution
        [0.05],   # transition: wider
        [0.01],   # subtelomeric: tight around zero
    ])

    return model


def label_windows(
    density_vector: np.ndarray,
    hmm: Optional["GaussianHMM"] = None,
    window_size: int = 6,
) -> Optional[np.ndarray]:
    """Run Viterbi decoding on windowed motif density.

    Args:
        density_vector: Per-base motif density from compute_motif_density.
        hmm: Trained or initialized GaussianHMM.
        window_size: Window size for downsampling density to HMM observations.

    Returns:
        State labels per window (0=telo, 1=transition, 2=subtelo),
        or None if HMM not available.
    """
    if hmm is None or not HAS_HMMLEARN:
        return None

    if len(density_vector) < window_size * 3:
        return None

    # Downsample density to window-level observations
    n_windows = len(density_vector) // window_size
    observations = np.array([
        density_vector[i * window_size:(i + 1) * window_size].mean()
        for i in range(n_windows)
    ]).reshape(-1, 1)

    if len(observations) < 3:
        return None

    try:
        _, state_path = hmm.decode(observations, algorithm="viterbi")
        return np.array(state_path)
    except Exception as e:
        logger.debug("HMM decoding failed: %s", e)
        return None


def hmm_boundary(
    state_path: Optional[np.ndarray],
    window_size: int = 6,
) -> Optional[int]:
    """Find the telomere-subtelomere boundary from HMM state path.

    Looks for the last transition from telomeric (state 0) to
    non-telomeric (state 1 or 2).

    Returns:
        Boundary position in base-pair coordinates, or None.
    """
    if state_path is None or len(state_path) == 0:
        return None

    # Find the last position where state is 0 (telomeric)
    telo_positions = np.where(state_path == 0)[0]
    if len(telo_positions) == 0:
        return None

    last_telo_window = telo_positions[-1]

    # If the entire path is telomeric, no boundary
    if last_telo_window >= len(state_path) - 1:
        return None

    # Convert window index to base position
    boundary_bp = (last_telo_window + 1) * window_size
    return boundary_bp


def expand_state_path(
    state_path: np.ndarray,
    window_size: int,
    target_length: int,
) -> np.ndarray:
    """Expand window-level state path to per-base resolution."""
    expanded = np.repeat(state_path, window_size)
    # Trim or pad to match target length
    if len(expanded) > target_length:
        expanded = expanded[:target_length]
    elif len(expanded) < target_length:
        expanded = np.pad(expanded, (0, target_length - len(expanded)),
                         constant_values=expanded[-1])
    return expanded


def fit_hmm(
    density_vectors: list[np.ndarray],
    window_size: int = 6,
    n_states: int = 3,
    n_iter: int = 100,
    random_seed: int = 42,
) -> Optional["GaussianHMM"]:
    """Train an HMM on a batch of motif-density vectors.

    Use this instead of :func:`build_hmm` when you have enough reads to
    reliably estimate emission and transition parameters from data.  The
    model is initialised with the same biologically-informed priors as
    ``build_hmm`` but all parameters (start, transition, means, covariances)
    are then refined via Baum-Welch EM on the supplied density vectors.

    Args:
        density_vectors: List of per-base motif density arrays, one per read.
        window_size: Window size for down-sampling density to HMM observations.
        n_states: Number of HMM states (default 3).
        n_iter: Maximum EM iterations.
        random_seed: For reproducibility.

    Returns:
        A trained GaussianHMM, or None if hmmlearn is not available or
        training fails.
    """
    if not HAS_HMMLEARN:
        return None

    # Start from biologically-informed priors
    model = build_hmm(n_states=n_states, random_seed=random_seed)
    if model is None:
        return None

    # Enable training of all parameter groups
    model.params = "stmc"
    model.n_iter = n_iter

    # Build observation sequences
    sequences = []
    lengths = []
    for dv in density_vectors:
        n_windows = len(dv) // window_size
        if n_windows < 3:
            continue
        obs = np.array([
            dv[i * window_size:(i + 1) * window_size].mean()
            for i in range(n_windows)
        ]).reshape(-1, 1)
        sequences.append(obs)
        lengths.append(len(obs))

    if not sequences:
        logger.warning("No density vectors long enough for HMM training")
        return None

    X = np.concatenate(sequences, axis=0)

    try:
        model.fit(X, lengths)
        logger.info(
            "HMM trained on %d sequences (%d total observations)",
            len(sequences), len(X),
        )
        return model
    except Exception as e:
        logger.error("HMM training failed: %s", e)
        return None
