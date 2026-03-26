"""Tests for chromosome arm assignment."""

import pytest

from teloprobe.anchor.arm_assigner import (
    parse_reference_name,
    _normalize_haplotype,
)


class TestParseReferenceName:
    def test_t2t_paternal(self):
        chrom, hap, arm = parse_reference_name("chr1_PATERNAL_P")
        assert chrom == "chr1"
        assert hap == "pat"
        assert arm == "p"

    def test_t2t_maternal(self):
        chrom, hap, arm = parse_reference_name("chr22_maternal_q")
        assert chrom == "chr22"
        assert hap == "mat"
        assert arm == "q"

    def test_simple_arm(self):
        chrom, hap, arm = parse_reference_name("chr1p")
        assert chrom == "chr1"
        assert arm == "p"
        assert hap is None

    def test_chrx(self):
        chrom, hap, arm = parse_reference_name("chrX_pat_q")
        assert chrom == "chrx"
        assert hap == "pat"
        assert arm == "q"

    def test_unparseable(self):
        chrom, hap, arm = parse_reference_name("random_contig")
        assert chrom is None
        assert hap is None
        assert arm is None

    def test_hap1_hap2(self):
        chrom, hap, arm = parse_reference_name("chr5_hap1_p")
        assert chrom == "chr5"
        assert hap == "pat"
        assert arm == "p"

    def test_underscore_arm(self):
        chrom, hap, arm = parse_reference_name("chr13_q")
        assert chrom == "chr13"
        assert arm == "q"


class TestNormalizeHaplotype:
    def test_paternal(self):
        assert _normalize_haplotype("PATERNAL") == "pat"
        assert _normalize_haplotype("pat") == "pat"
        assert _normalize_haplotype("hap1") == "pat"

    def test_maternal(self):
        assert _normalize_haplotype("MATERNAL") == "mat"
        assert _normalize_haplotype("mat") == "mat"
        assert _normalize_haplotype("hap2") == "mat"

    def test_unknown(self):
        assert _normalize_haplotype("other") == "unknown"
