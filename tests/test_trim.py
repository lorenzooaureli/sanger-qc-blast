"""Tests for trimming algorithms."""

import pytest
from sanger_qc_trim.trim import trim_mott, trim_ends, apply_trim


class TestTrimMott:
    """Tests for Mott/Kadane trimming algorithm."""

    def test_basic_trimming(self):
        """Test basic Mott trimming on simple quality array."""
        # quals = [10, 10, 30, 30, 30, 10], T=20
        # weights = [-10, -10, 10, 10, 10, -10]
        # Best window: indices 2-5 (length 3)
        quals = [10, 10, 30, 30, 30, 10]
        result = trim_mott(quals, 20)
        assert result == (2, 5)

    def test_all_high_quality(self):
        """Test when all bases are high quality."""
        quals = [30, 30, 30, 30]
        result = trim_mott(quals, 20)
        assert result == (0, 4)

    def test_all_low_quality(self):
        """Test when all bases are low quality."""
        quals = [10, 10, 10, 10]
        result = trim_mott(quals, 20)
        assert result == (0, 0)

    def test_empty_input(self):
        """Test with empty quality list."""
        quals = []
        result = trim_mott(quals, 20)
        assert result == (0, 0)

    def test_single_high_quality(self):
        """Test with single high quality base."""
        quals = [30]
        result = trim_mott(quals, 20)
        assert result == (0, 1)

    def test_low_then_high(self):
        """Test low quality followed by high quality."""
        quals = [10, 10, 30, 30]
        result = trim_mott(quals, 20)
        assert result == (2, 4)

    def test_high_then_low(self):
        """Test high quality followed by low quality."""
        quals = [30, 30, 10, 10]
        result = trim_mott(quals, 20)
        assert result == (0, 2)

    def test_mixed_quality(self):
        """Test mixed quality with multiple regions."""
        # Two high quality regions: prefer the longer/higher scoring one
        quals = [30, 30, 10, 10, 35, 35, 35]
        result = trim_mott(quals, 20)
        # Second region (indices 4-7) has better score
        assert result == (4, 7)


class TestTrimEnds:
    """Tests for ends trimming algorithm."""

    def test_basic_trimming(self):
        """Test basic ends trimming."""
        quals = [15, 25, 25, 15]
        result = trim_ends(quals, 20)
        assert result == (1, 3)

    def test_all_high_quality(self):
        """Test when all bases are high quality."""
        quals = [30, 30, 30, 30]
        result = trim_ends(quals, 20)
        assert result == (0, 4)

    def test_all_low_quality(self):
        """Test when all bases are low quality."""
        quals = [10, 10, 10, 10]
        result = trim_ends(quals, 20)
        assert result == (0, 0)

    def test_empty_input(self):
        """Test with empty quality list."""
        quals = []
        result = trim_ends(quals, 20)
        assert result == (0, 0)

    def test_single_high_quality(self):
        """Test with single high quality base."""
        quals = [30]
        result = trim_ends(quals, 20)
        assert result == (0, 1)

    def test_first_base_only(self):
        """Test when only first base passes threshold."""
        quals = [30, 10, 10, 10]
        result = trim_ends(quals, 20)
        assert result == (0, 1)

    def test_last_base_only(self):
        """Test when only last base passes threshold."""
        quals = [10, 10, 10, 30]
        result = trim_ends(quals, 20)
        assert result == (3, 4)

    def test_middle_bases_low(self):
        """Test when middle bases are low but ends are high."""
        quals = [30, 10, 10, 30]
        result = trim_ends(quals, 20)
        # Keeps from first high to last high
        assert result == (0, 4)


class TestApplyTrim:
    """Tests for apply_trim function."""

    def test_apply_mott(self):
        """Test applying Mott trimming to sequence."""
        seq = "ACGTACGT"
        quals = [10, 10, 30, 30, 30, 10, 10, 10]
        trimmed_seq, trimmed_quals, start, end = apply_trim(seq, quals, "mott", 20)

        assert start == 2
        assert end == 5
        assert trimmed_seq == "GTA"
        assert trimmed_quals == [30, 30, 30]

    def test_apply_ends(self):
        """Test applying ends trimming to sequence."""
        seq = "ACGT"
        quals = [15, 25, 25, 15]
        trimmed_seq, trimmed_quals, start, end = apply_trim(seq, quals, "ends", 20)

        assert start == 1
        assert end == 3
        assert trimmed_seq == "CG"
        assert trimmed_quals == [25, 25]

    def test_invalid_method(self):
        """Test with invalid trimming method."""
        seq = "ACGT"
        quals = [20, 20, 20, 20]

        with pytest.raises(ValueError, match="Unknown trimming method"):
            apply_trim(seq, quals, "invalid", 20)

    def test_empty_result(self):
        """Test when trimming results in empty sequence."""
        seq = "ACGT"
        quals = [10, 10, 10, 10]
        trimmed_seq, trimmed_quals, start, end = apply_trim(seq, quals, "mott", 20)

        assert start == 0
        assert end == 0
        assert trimmed_seq == ""
        assert trimmed_quals == []


class TestQualityThresholds:
    """Tests for different quality thresholds."""

    def test_q20_threshold(self):
        """Test with Q20 threshold."""
        quals = [20, 20, 10, 30]
        result = trim_mott(quals, 20)
        # [0, 0, -10, 10] -> best is index 3-4
        assert result == (3, 4)

    def test_q30_threshold(self):
        """Test with Q30 threshold."""
        quals = [20, 20, 35, 35]
        result = trim_mott(quals, 30)
        # [-10, -10, 5, 5] -> best is indices 2-4
        assert result == (2, 4)

    def test_boundary_quality(self):
        """Test with quality exactly at threshold."""
        quals = [20, 20, 20]
        result = trim_mott(quals, 20)
        # [0, 0, 0] -> no positive weights, returns (0, 0)
        assert result == (0, 0)
