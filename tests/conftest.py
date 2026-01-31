"""Pytest configuration and fixtures."""

import pytest
from pathlib import Path


@pytest.fixture
def tmp_output_dir(tmp_path):
    """Create a temporary output directory."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return output_dir


@pytest.fixture
def sample_sequence():
    """Sample DNA sequence for testing."""
    return "ACGTACGTACGTACGT"


@pytest.fixture
def sample_qualities():
    """Sample quality scores for testing."""
    return [20, 25, 30, 35, 30, 25, 20, 15, 10, 15, 20, 25, 30, 35, 30, 25]
