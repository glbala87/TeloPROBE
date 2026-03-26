"""Tests for TVR decomposition."""

import pytest

from teloprobe.candidate.tvr_decomposer import (
    decompose_tvr,
    tvr_fraction,
    tvr_profile_vector,
    tvr_pattern_string,
)


class TestDecomposeTVR:
    def test_canonical_only(self):
        seq = "TAACCC" * 100
        comp = decompose_tvr(seq, boundary=600)
        assert "TTAGGG" in comp  # canonical mapped to TTAGGG
        assert comp["TTAGGG"] == 100

    def test_with_variants(self):
        seq = "TAACCC" * 80 + "TGAGGG" * 20
        comp = decompose_tvr(seq)
        assert comp.get("TTAGGG", 0) > 0
        assert comp.get("TGAGGG", 0) > 0

    def test_boundary_limits_region(self):
        seq = "TAACCC" * 50 + "ATCGATCG" * 50
        comp_full = decompose_tvr(seq)
        comp_bounded = decompose_tvr(seq, boundary=300)
        assert comp_bounded["TTAGGG"] == 50

    def test_empty_sequence(self):
        comp = decompose_tvr("")
        assert len(comp) == 0


class TestTVRFraction:
    def test_all_canonical(self):
        comp = {"TTAGGG": 100}
        assert tvr_fraction(comp) == 0.0

    def test_mixed(self):
        comp = {"TTAGGG": 80, "TGAGGG": 20}
        assert tvr_fraction(comp) == pytest.approx(0.2)

    def test_empty(self):
        assert tvr_fraction({}) == 0.0


class TestTVRPatternString:
    def test_output_format(self):
        comp = {"TTAGGG": 100, "TGAGGG": 5}
        s = tvr_pattern_string(comp)
        assert "TTAGGG:100" in s
        assert "TGAGGG:5" in s

    def test_empty(self):
        assert tvr_pattern_string({}) == "none"


class TestTVRProfileVector:
    def test_length(self):
        comp = {"TTAGGG": 100}
        from teloprobe.constants import VARIANT_MOTIFS
        vec = tvr_profile_vector(comp)
        assert len(vec) == len(VARIANT_MOTIFS)

    def test_normalized(self):
        comp = {"TTAGGG": 100, "TGAGGG": 50}
        vec = tvr_profile_vector(comp)
        # All values should be between 0 and 1
        assert all(0 <= v <= 1 for v in vec)
