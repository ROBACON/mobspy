"""Shared fixtures and utilities for MobsPy tests."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

# Root of the repository
ROOT = Path(__file__).resolve().parent.parent
TEST_TOOLS = ROOT / "test_tools"


def compare_model(comp_results: str, file_name: str) -> bool:
    """Compare compiled model output against an expected output file.

    Normalises whitespace before comparing line-by-line.
    """
    path = (
        TEST_TOOLS / file_name if not Path(file_name).is_absolute() else Path(file_name)
    )
    expected_lines = path.read_text(encoding="utf-8").splitlines()
    result_lines = comp_results.splitlines()

    def _normalise(text: str) -> str:
        return re.sub(r"\s+", " ", text.strip())

    for result_line, expected_line in zip(result_lines, expected_lines, strict=False):
        if _normalise(result_line) != _normalise(expected_line):
            return False
    return True


def compare_model_ignore_order(comp_results: str, file_name: str) -> bool:
    """Compare compiled model output ignoring line order."""
    path = (
        TEST_TOOLS / file_name if not Path(file_name).is_absolute() else Path(file_name)
    )
    expected_lines = {
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    result_lines = {line.strip() for line in comp_results.splitlines() if line.strip()}
    return expected_lines == result_lines


@pytest.fixture()
def test_tools_dir() -> Path:
    """Return the path to the test_tools directory."""
    return TEST_TOOLS
