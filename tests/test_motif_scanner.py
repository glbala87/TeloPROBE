"""Tests for motif scanning module."""

import numpy as np
import pytest

from teloprobe.candidate.motif_scanner import (
    scan_motifs,
    compute_motif_density,
    compute_motif_density_fast,
    quick_filter,
    longest_telomeric_run,
    telomeric_fraction,
    entropy_score,
)
from tests.conftest import make_telomeric_sequence


class TestScanMotifs:
    def test_canonical_motif_detection(self):
        """TAACCC repeats should be fully detected."""
        seq = "TAACCC" * 50
        hits = scan_motifs(seq)
        assert hits.sum() == len(seq)  # every position covered

    def test_no_motifs_in_random_sequence(self):
        """Random sequence should have very few/no motif hits."""
        seq = "ATCGATCGATCGATCG" * 20
        hits = scan_motifs(seq)
        assert hits.sum() / len(seq) < 0.1

    def test_mixed_sequence(self):
        """Telomeric region should be detected in mixed sequence."""
        seq = make_telomeric_sequence(600, 400)
        hits = scan_motifs(seq)

        # First 600bp should be mostly hits
        telo_density = hits[:600].sum() / 600
        assert telo_density > 0.8

        # Last 400bp should be mostly non-hits
        subtelo_density = hits[600:].sum() / 400
        assert subtelo_density < 0.3

    def test_empty_sequence(self):
        hits = scan_motifs("")
        assert len(hits) == 0

    def test_short_sequence(self):
        hits = scan_motifs("ACGT")
        assert len(hits) == 4
        assert hits.sum() == 0

    def test_variant_motifs(self):
        """TVR motifs should be detected when include_variants=True."""
        seq = "TGAGGG" * 20  # variant repeat
        hits_with = scan_motifs(seq, include_variants=True)
        hits_without = scan_motifs(seq, include_variants=False)
        assert hits_with.sum() > hits_without.sum()


class TestMotifDensity:
    def test_pure_telomere_density(self):
        seq = "TAACCC" * 100
        hits = scan_motifs(seq)
        density = compute_motif_density(hits, window_size=60)
        assert np.mean(density) > 0.9

    def test_density_transition(self):
        """Density should drop at telomere-subtelomere boundary."""
        seq = make_telomeric_sequence(500, 500)
        hits = scan_motifs(seq)
        density = compute_motif_density(hits, window_size=60)

        before = np.mean(density[100:400])
        after = np.mean(density[600:900])
        assert before > 0.7
        assert after < 0.3

    def test_fast_density_matches_regular(self):
        """Fast density should be close to regular density."""
        seq = make_telomeric_sequence(300, 300)
        hits = scan_motifs(seq)
        d1 = compute_motif_density(hits, window_size=60)
        d2 = compute_motif_density_fast(hits, window_size=60)
        np.testing.assert_allclose(d1, d2, atol=0.15)


class TestQuickFilter:
    def test_telomeric_passes(self):
        seq = "TAACCC" * 200
        assert quick_filter(seq, min_repeats=100) is True

    def test_random_fails(self):
        seq = "ATCGATCG" * 200
        assert quick_filter(seq, min_repeats=100) is False

    def test_borderline(self):
        seq = "TAACCC" * 50  # exactly 50 repeats
        assert quick_filter(seq, min_repeats=50) is True
        assert quick_filter(seq, min_repeats=51) is False


class TestHelpers:
    def test_longest_run(self):
        arr = np.array([1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0])
        assert longest_telomeric_run(arr) == 5

    def test_longest_run_empty(self):
        assert longest_telomeric_run(np.array([])) == 0

    def test_telomeric_fraction(self):
        arr = np.array([1, 1, 1, 0, 0])
        assert telomeric_fraction(arr) == pytest.approx(0.6)

    def test_entropy_low_complexity(self):
        seq = "AAAAAA" * 20
        e = entropy_score(seq)
        assert e < 0.5

    def test_entropy_high_complexity(self):
        seq = "ACGT" * 30
        e = entropy_score(seq)
        assert e > 1.5
