"""Tests for confidence scoring and bootstrap CIs."""

import numpy as np
import pytest

from teloprobe.estimate.confidence import (
    bootstrap_ci,
    compute_confidence_score,
    evidence_weights,
)
from teloprobe.constants import Mode
from teloprobe.config import Config
from teloprobe.models.telomere_call import TelomereCall, EvidenceLayers


class TestBootstrapCI:
    def test_basic_ci(self):
        rng = np.random.RandomState(42)
        values = rng.normal(1000, 100, 200)
        ci_low, ci_high = bootstrap_ci(values, n_boot=500)
        assert ci_low < 1000
        assert ci_high > 1000
        assert ci_low < ci_high

    def test_narrow_ci_with_tight_data(self):
        values = np.full(100, 500.0) + np.random.RandomState(42).normal(0, 1, 100)
        ci_low, ci_high = bootstrap_ci(values, n_boot=500)
        assert ci_high - ci_low < 10

    def test_single_value(self):
        values = np.array([42.0])
        ci_low, ci_high = bootstrap_ci(values)
        assert ci_low == pytest.approx(42.0)

    def test_two_values(self):
        values = np.array([40.0, 60.0])
        ci_low, ci_high = bootstrap_ci(values)
        assert ci_low <= 50
        assert ci_high >= 50


class TestConfidenceScore:
    def test_high_confidence(self):
        call = TelomereCall(
            read_id="test", read_length=1000,
            evidence=EvidenceLayers(motif_score=0.95, boundary_score=1.0, anchor_score=0.9),
        )
        config = Config()
        config.mode = Mode.TELOSEQ
        score = compute_confidence_score(call, config)
        assert score > 0.8

    def test_low_confidence(self):
        call = TelomereCall(
            read_id="test", read_length=1000,
            evidence=EvidenceLayers(motif_score=0.3, boundary_score=0.2, anchor_score=0.0),
        )
        config = Config()
        score = compute_confidence_score(call, config)
        assert score < 0.5

    def test_failed_qc_penalized(self):
        from teloprobe.constants import QCStatus
        call = TelomereCall(
            read_id="test", read_length=1000,
            qc_status=QCStatus.TOO_FEW_REPEATS,
            evidence=EvidenceLayers(motif_score=0.9, boundary_score=0.9, anchor_score=0.9),
        )
        config = Config()
        score = compute_confidence_score(call, config)
        assert score < 0.2


class TestEvidenceWeights:
    def test_teloprobe_mode(self):
        w = evidence_weights(Mode.TELOSEQ)
        assert w["motif"] == pytest.approx(0.4)
        assert w["boundary"] == pytest.approx(0.4)
        assert w["anchor"] == pytest.approx(0.2)

    def test_wgs_mode(self):
        w = evidence_weights(Mode.WGS)
        assert w["anchor"] > w["motif"]  # anchor more important in WGS

    def test_weights_sum_to_one(self):
        for mode in Mode:
            w = evidence_weights(mode)
            assert sum(w.values()) == pytest.approx(1.0)
