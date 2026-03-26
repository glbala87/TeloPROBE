"""Tests for HMM labeler module."""

import numpy as np
import pytest

from teloprobe.segment.hmm_labeler import (
    build_hmm,
    label_windows,
    hmm_boundary,
    expand_state_path,
)


class TestBuildHMM:
    def test_creates_hmm(self):
        hmm = build_hmm()
        assert hmm is not None
        assert hmm.n_components == 3

    def test_start_probabilities(self):
        hmm = build_hmm()
        assert hmm.startprob_[0] > 0.8  # mostly starts in telo state

    def test_transition_matrix_valid(self):
        hmm = build_hmm()
        # Rows should sum to 1
        for row in hmm.transmat_:
            assert abs(sum(row) - 1.0) < 1e-6

    def test_custom_density(self):
        hmm = build_hmm(telo_mean=0.9, subtelo_mean=0.1)
        assert hmm.means_[0][0] == pytest.approx(0.9)
        assert hmm.means_[2][0] == pytest.approx(0.1)


class TestLabelWindows:
    def test_labels_telomeric_region(self):
        hmm = build_hmm()
        # High density followed by low density
        density = np.concatenate([
            np.full(300, 0.9),
            np.full(300, 0.05),
        ])
        labels = label_windows(density, hmm, window_size=6)
        assert labels is not None

        # First windows should be state 0 (telomeric)
        assert labels[0] == 0
        # Last windows should be state 2 (subtelomeric)
        assert labels[-1] == 2

    def test_short_signal_returns_none(self):
        hmm = build_hmm()
        density = np.array([0.9, 0.1])
        labels = label_windows(density, hmm, window_size=6)
        assert labels is None

    def test_none_hmm_returns_none(self):
        density = np.full(600, 0.5)
        labels = label_windows(density, None)
        assert labels is None


class TestHMMBoundary:
    def test_finds_boundary(self):
        # State path: telo telo telo telo subtelo subtelo subtelo subtelo
        path = np.array([0, 0, 0, 0, 2, 2, 2, 2])
        bnd = hmm_boundary(path, window_size=6)
        assert bnd is not None
        # Last telo window is index 3, boundary = (3+1)*6 = 24
        assert bnd == 24

    def test_no_boundary_all_telo(self):
        path = np.array([0, 0, 0, 0, 0])
        bnd = hmm_boundary(path, window_size=6)
        assert bnd is None

    def test_none_path(self):
        assert hmm_boundary(None) is None

    def test_empty_path(self):
        assert hmm_boundary(np.array([])) is None


class TestExpandStatePath:
    def test_expansion(self):
        path = np.array([0, 0, 2, 2])
        expanded = expand_state_path(path, window_size=3, target_length=12)
        assert len(expanded) == 12
        assert expanded[0] == 0
        assert expanded[-1] == 2

    def test_padding(self):
        path = np.array([0, 2])
        expanded = expand_state_path(path, window_size=3, target_length=10)
        assert len(expanded) == 10
