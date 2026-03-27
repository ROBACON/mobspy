"""Tests for multiline reaction syntax."""

from __future__ import annotations

from mobspy import BaseSpecies, Simulation


def test_multiline():
    A, B = BaseSpecies(2)
    # fmt: off
    (
        A >>
        B [1]
    )
    # fmt: on
    MySim = Simulation(A | B)
    MySim.level = -1
    result = MySim.compile()
    assert result is not None
    assert len(result) > 0


def test_ruff_style_multiline_reaction():
    A, B, C, D = BaseSpecies(4)
    # fmt: off
    (
        A + B >>
        C + D [1]
    )
    # fmt: on
    MySim = Simulation(A | B | C | D)
    MySim.level = -1
    result = MySim.compile()
    assert result is not None
    assert len(result) > 0


def test_simple_multiline_reaction():
    A, B = BaseSpecies(2)
    # fmt: off
    (
        A >>
        B [1]
    )
    # fmt: on
    MySim = Simulation(A | B)
    MySim.level = -1
    result = MySim.compile()
    assert result is not None
    assert len(result) > 0


def test_complex_multiline_reaction():
    A, B, C, D = BaseSpecies(4)
    # fmt: off
    (
        2*A + B >>
        3*C + 2*D [1]
    )
    # fmt: on
    MySim = Simulation(A | B | C | D)
    MySim.level = -1
    result = MySim.compile()
    assert result is not None
    assert len(result) > 0


def test_normal_reaction():
    A, B, C = BaseSpecies(3)
    A + B >> C[1]
    MySim = Simulation(A | B | C)
    MySim.level = -1
    result = MySim.compile()
    assert result is not None
    assert len(result) > 0
