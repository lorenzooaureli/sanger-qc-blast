"""Tests for QC statistics computation."""

import pytest
from sanger_qc_trim.qc import compute_qc_metrics, compute_summary_stats, _longest_hq_stretch


class TestQCMetrics:
    """Tests for QC metrics computation."""

    def test_basic_metrics(self):
        """Test basic QC metrics computation."""
        metrics = compute_qc_metrics(
            sample_id="test_sample",
            source_file="/path/to/test.ab1",
            file_format="ab1",
            seq="ACGTACGT",
            quals=[20, 20, 30, 30, 20, 20, 10, 10],
            trim_start=0,
            trim_end=6,
            qthreshold=20,
            min_length=5,
        )

        assert metrics["sample_id"] == "test_sample"
        assert metrics["source_file"] == "/path/to/test.ab1"
        assert metrics["format"] == "ab1"
        assert metrics["raw_length"] == 8
        assert metrics["trimmed_length"] == 6
        assert metrics["passed_minlen"] == "yes"

    def test_mean_quality(self):
        """Test mean quality calculation."""
        metrics = compute_qc_metrics(
            sample_id="test",
            source_file="test.ab1",
            file_format="ab1",
            seq="ACGT",
            quals=[10, 20, 30, 40],
            trim_start=0,
            trim_end=4,
            qthreshold=20,
            min_length=1,
        )

        # Mean = (10 + 20 + 30 + 40) / 4 = 25
        assert metrics["mean_q"] == 25.0
        assert metrics["median_q"] == 25.0

    def test_q20_q30_percentages(self):
        """Test Q20 and Q30 percentage calculations."""
        # quals = [20, 20, 10, 30] -> 3/4 are Q20+, 1/4 are Q30+
        metrics = compute_qc_metrics(
            sample_id="test",
            source_file="test.ab1",
            file_format="ab1",
            seq="ACGT",
            quals=[20, 20, 10, 30],
            trim_start=0,
            trim_end=4,
            qthreshold=20,
            min_length=1,
        )

        assert metrics["pct_q20"] == 0.75
        assert metrics["pct_q30"] == 0.25

    def test_gc_content(self):
        """Test GC content calculation."""
        # seq = "ACGT" -> 2 GC out of 4 = 50%
        metrics = compute_qc_metrics(
            sample_id="test",
            source_file="test.ab1",
            file_format="ab1",
            seq="ACGT",
            quals=[20, 20, 20, 20],
            trim_start=0,
            trim_end=4,
            qthreshold=20,
            min_length=1,
        )

        assert metrics["gc_percent"] == 50.0

    def test_gc_content_with_n(self):
        """Test GC content with N bases."""
        # seq = "ACGTN" -> 2 GC out of 4 non-N = 50%
        metrics = compute_qc_metrics(
            sample_id="test",
            source_file="test.ab1",
            file_format="ab1",
            seq="ACGTN",
            quals=[20, 20, 20, 20, 20],
            trim_start=0,
            trim_end=5,
            qthreshold=20,
            min_length=1,
        )

        assert metrics["gc_percent"] == 50.0
        assert metrics["n_count"] == 1

    def test_expected_errors(self):
        """Test expected errors calculation."""
        # For Q20: 10^(-20/10) = 0.01
        # For Q30: 10^(-30/10) = 0.001
        metrics = compute_qc_metrics(
            sample_id="test",
            source_file="test.ab1",
            file_format="ab1",
            seq="AC",
            quals=[20, 30],
            trim_start=0,
            trim_end=2,
            qthreshold=20,
            min_length=1,
        )

        # Expected: 0.01 + 0.001 = 0.011
        assert abs(metrics["expected_errors"] - 0.011) < 0.001

    def test_passed_minlen_yes(self):
        """Test passed_minlen flag when passing."""
        metrics = compute_qc_metrics(
            sample_id="test",
            source_file="test.ab1",
            file_format="ab1",
            seq="ACGTACGT",
            quals=[20] * 8,
            trim_start=0,
            trim_end=8,
            qthreshold=20,
            min_length=5,
        )

        assert metrics["passed_minlen"] == "yes"

    def test_passed_minlen_no(self):
        """Test passed_minlen flag when failing."""
        metrics = compute_qc_metrics(
            sample_id="test",
            source_file="test.ab1",
            file_format="ab1",
            seq="ACGT",
            quals=[20] * 4,
            trim_start=0,
            trim_end=4,
            qthreshold=20,
            min_length=10,
        )

        assert metrics["passed_minlen"] == "no"


class TestLongestHQStretch:
    """Tests for longest high-quality stretch calculation."""

    def test_single_stretch(self):
        """Test single continuous high-quality stretch."""
        quals = [30, 30, 30, 10]
        result = _longest_hq_stretch(quals, 20)
        assert result == 3

    def test_multiple_stretches(self):
        """Test multiple high-quality stretches."""
        quals = [30, 30, 10, 30, 30, 30, 30]
        result = _longest_hq_stretch(quals, 20)
        assert result == 4

    def test_all_high_quality(self):
        """Test when all bases are high quality."""
        quals = [30, 30, 30, 30]
        result = _longest_hq_stretch(quals, 20)
        assert result == 4

    def test_no_high_quality(self):
        """Test when no bases are high quality."""
        quals = [10, 10, 10, 10]
        result = _longest_hq_stretch(quals, 20)
        assert result == 0

    def test_empty_input(self):
        """Test with empty input."""
        quals = []
        result = _longest_hq_stretch(quals, 20)
        assert result == 0


class TestSummaryStats:
    """Tests for summary statistics computation."""

    def test_basic_summary(self):
        """Test basic summary statistics."""
        metrics_list = [
            {
                "raw_length": 100,
                "trimmed_length": 80,
                "mean_q": 25.0,
                "pct_q20": 0.8,
                "pct_q30": 0.5,
                "gc_percent": 50.0,
                "expected_errors": 1.0,
                "passed_minlen": "yes",
            },
            {
                "raw_length": 120,
                "trimmed_length": 90,
                "mean_q": 30.0,
                "pct_q20": 0.9,
                "pct_q30": 0.7,
                "gc_percent": 45.0,
                "expected_errors": 0.5,
                "passed_minlen": "yes",
            },
        ]

        summary = compute_summary_stats(metrics_list)

        assert summary["total_reads"] == 2
        assert summary["mean_raw_length"] == 110.0
        assert summary["mean_trimmed_length"] == 85.0
        assert summary["mean_mean_q"] == 27.5
        assert summary["reads_passed_minlen"] == 2
        assert summary["reads_failed_minlen"] == 0
        assert summary["pct_passed"] == 100.0

    def test_summary_with_failures(self):
        """Test summary with some failed reads."""
        metrics_list = [
            {
                "raw_length": 100,
                "trimmed_length": 80,
                "mean_q": 25.0,
                "pct_q20": 0.8,
                "pct_q30": 0.5,
                "gc_percent": 50.0,
                "expected_errors": 1.0,
                "passed_minlen": "yes",
            },
            {
                "raw_length": 100,
                "trimmed_length": 20,
                "mean_q": 15.0,
                "pct_q20": 0.3,
                "pct_q30": 0.1,
                "gc_percent": 40.0,
                "expected_errors": 5.0,
                "passed_minlen": "no",
            },
        ]

        summary = compute_summary_stats(metrics_list)

        assert summary["total_reads"] == 2
        assert summary["reads_passed_minlen"] == 1
        assert summary["reads_failed_minlen"] == 1
        assert summary["pct_passed"] == 50.0

    def test_empty_metrics_list(self):
        """Test summary with empty metrics list."""
        summary = compute_summary_stats([])

        assert summary["total_reads"] == 0
        assert summary["mean_raw_length"] == 0.0
        assert summary["reads_passed_minlen"] == 0
        assert summary["reads_failed_minlen"] == 0
        assert summary["pct_passed"] == 0.0
