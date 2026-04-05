"""Tests proving models in separate functions work correctly.

Each function defines a complete model, compiles it, and runs it.
The results must be correct regardless of what other functions did before.
"""

from __future__ import annotations

from mobspy import BaseSpecies, New, Simulation, Zero


def _simple_decay() -> Simulation:
    A = BaseSpecies(["A"])
    A >> Zero @ 1.0
    A(100)
    S = Simulation(A)
    S.duration = 5
    S.level = -1
    S.plot_data = False
    return S


def _bimolecular() -> Simulation:
    A, B = BaseSpecies(["A", "B"])
    A + B >> Zero @ 0.01
    A(10), B(5)
    S = Simulation(A | B)
    S.method = "stochastic"
    S.duration = (A <= 0) | (B <= 0)
    S.level = -1
    S.plot_data = False
    return S


def _inherited_decay() -> Simulation:
    Parent = BaseSpecies(["Parent"])
    C1, C2 = New(Parent, 2)
    Parent >> Zero @ 0.5
    C1(80)
    C2(40)
    S = Simulation(C1 | C2)
    S.duration = 5
    S.level = -1
    S.plot_data = False
    return S


def _reversible() -> Simulation:
    A, B = BaseSpecies(["A", "B"])
    A >> B @ (1.0, 0.5)
    A(100)
    S = Simulation(A | B)
    S.duration = 10
    S.level = -1
    S.plot_data = False
    return S


class TestFunctionIsolation:
    def test_decay_then_bimolecular(self) -> None:
        S1 = _simple_decay()
        S1.run(plot_data=False)
        assert S1.fres["A"][-1] < 1

        S2 = _bimolecular()
        S2.run(plot_data=False)
        assert S2.fres["A"][-1] == 0 or S2.fres["B"][-1] == 0

    def test_bimolecular_then_decay(self) -> None:
        S1 = _bimolecular()
        S1.run(plot_data=False)
        assert S1.fres["A"][-1] == 0 or S1.fres["B"][-1] == 0

        S2 = _simple_decay()
        S2.run(plot_data=False)
        assert S2.fres["A"][-1] < 1

    def test_inherited_then_simple(self) -> None:
        S1 = _inherited_decay()
        S1.run(plot_data=False)
        assert S1.fres["C1"][-1] < 10
        assert S1.fres["C2"][-1] < 10

        S2 = _simple_decay()
        S2.run(plot_data=False)
        assert S2.fres["A"][-1] < 1

    def test_reversible_then_inherited(self) -> None:
        S1 = _reversible()
        S1.run(plot_data=False)
        a_final = S1.fres["A"][-1]
        b_final = S1.fres["B"][-1]
        assert a_final + b_final > 99
        assert b_final > 10

        S2 = _inherited_decay()
        S2.run(plot_data=False)
        assert S2.fres["C1"][-1] < 10

    def test_five_sequential_models(self) -> None:
        for _ in range(5):
            S = _simple_decay()
            S.run(plot_data=False)
            assert S.fres["A"][-1] < 1

    def test_same_names_different_rates(self) -> None:
        """Both functions define A >> Zero but with different rates.
        Each must get its own rate."""

        def fast():
            A = BaseSpecies(["A"])
            A >> Zero @ 10.0
            A(100)
            S = Simulation(A)
            S.duration = 1
            S.level = -1
            S.plot_data = False
            return S

        def slow():
            A = BaseSpecies(["A"])
            A >> Zero @ 0.01
            A(100)
            S = Simulation(A)
            S.duration = 1
            S.level = -1
            S.plot_data = False
            return S

        S_fast = fast()
        S_slow = slow()
        S_fast.run(plot_data=False)
        S_slow.run(plot_data=False)
        assert S_fast.fres["A"][-1] < 1
        assert S_slow.fres["A"][-1] > 90

    def test_same_names_different_counts(self) -> None:
        """Both define A with different initial counts."""

        def ten():
            A = BaseSpecies(["A"])
            A >> Zero @ 1.0
            A(10)
            S = Simulation(A)
            S.duration = 0.01
            S.level = -1
            S.plot_data = False
            return S

        def thousand():
            A = BaseSpecies(["A"])
            A >> Zero @ 1.0
            A(1000)
            S = Simulation(A)
            S.duration = 0.01
            S.level = -1
            S.plot_data = False
            return S

        S1 = ten()
        S2 = thousand()
        S1.run(plot_data=False)
        S2.run(plot_data=False)
        assert S1.fres["A"][0] == 10
        assert S2.fres["A"][0] == 1000

    def test_same_names_different_reactions(self) -> None:
        """One defines A >> B, the other defines A >> C.
        Each must only have its own reaction."""

        def to_b():
            A, B = BaseSpecies(["A", "B"])
            A >> B @ 1.0
            A(100)
            S = Simulation(A | B)
            S.duration = 5
            S.level = -1
            S.plot_data = False
            return S

        def to_c():
            A, C = BaseSpecies(["A", "C"])
            A >> C @ 1.0
            A(100)
            S = Simulation(A | C)
            S.duration = 5
            S.level = -1
            S.plot_data = False
            return S

        S1 = to_b()
        S2 = to_c()
        S1.run(plot_data=False)
        S2.run(plot_data=False)
        assert S1.fres["B"][-1] > 90
        assert S2.fres["C"][-1] > 90
