"""Tests for consensus boundary detection."""

import numpy as np
import pytest

from teloprobe.segment.boundary import (
    consensus_boundary,
    classify_boundary,
    determine_telomere_side,
)
from teloprobe.constants import QCStatus
from teloprobe.config import Config
from tests.conftest import make_telomeric_sequence


class TestConsensusBoundary:
    def test_both_agree(self):
        """When HMM and PELT agree, return consensus with high confidence."""
        motif = np.concatenate([np.ones(500), np.zeros(500)])
        bnd, method, score = consensus_boundary(500, 510, motif, tolerance=60)
        assert bnd is not None
        assert method == "consensus"
        assert score > 0.7
        assert abs(bnd - 505) <= 5

    def test_both_disagree(self):
        """When boundaries differ significantly, pick best contrast."""
        motif = np.concatenate([np.ones(500), np.zeros(500)])
        bnd, method, score = consensus_boundary(300, 500, motif, tolerance=60)
        assert bnd is not None
        assert score <= 0.6

    def test_hmm_only(self):
        motif = np.concatenate([np.ones(500), np.zeros(500)])
        bnd, method, score = consensus_boundary(500, None, motif)
        assert bnd == 500
        assert method == "hmm"
        assert score == pytest.approx(0.6)

    def test_pelt_only(self):
        motif = np.concatenate([np.ones(500), np.zeros(500)])
        bnd, method, score = consensus_boundary(None, 500, motif)
        assert bnd == 500
        assert method == "pelt"

    def test_neither(self):
        motif = np.zeros(1000)
        bnd, method, score = consensus_boundary(None, None, motif)
        # May find threshold boundary or None
        assert score <= 0.3


class TestClassifyBoundary:
    @pytest.fixture
    def config(self):
        c = Config()
        c.min_read_length_boundary = 50
        c.min_repeats = 10
        c.filter_width = 10
        c.start_window_frac = 0.3
        c.start_repeats_frac = 0.8
        c.min_qual_non_telo = 9
        return c

    def test_good_boundary(self, config):
        motif = np.concatenate([np.ones(500), np.zeros(500)])
        quals = np.full(1000, 20, dtype=np.int32)
        status = classify_boundary(500, 1000, motif, quals, 80, config)
        assert status == QCStatus.GOOD

    def test_too_short(self, config):
        motif = np.ones(30)
        status = classify_boundary(15, 30, motif, None, 5, config)
        assert status == QCStatus.TOO_SHORT

    def test_too_few_repeats(self, config):
        motif = np.concatenate([np.ones(30), np.zeros(500)])
        status = classify_boundary(30, 530, motif, None, 5, config)
        assert status == QCStatus.TOO_FEW_REPEATS

    def test_no_boundary(self, config):
        motif = np.ones(500)
        status = classify_boundary(None, 500, motif, None, 80, config)
        assert status == QCStatus.NO_BOUNDARY


class TestDetermineTeloSide:
    def test_5prime(self):
        motif = np.concatenate([np.ones(500), np.zeros(500)])
        side = determine_telomere_side(motif, 500, 1000)
        assert side == "5prime"

    def test_no_boundary_is_internal(self):
        motif = np.zeros(1000)
        side = determine_telomere_side(motif, None, 1000)
        assert side == "internal"
