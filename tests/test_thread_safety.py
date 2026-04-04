"""Tests verifying that MobsPy model definition is thread-safe.

Two threads define independent models concurrently and compile them.
Results must be correct and non-interfering.
"""

from __future__ import annotations

import threading
from typing import Any

from mobspy import BaseSpecies, New, Simulation, Zero


def _define_and_compile_model_a() -> dict[str, Any]:
    """Define a simple birth-death model and compile it."""
    A = BaseSpecies(["A"])
    _ = Zero >> A @ 1
    _ = A >> Zero @ 0.1
    A(100)
    sim = Simulation(A)
    sim.duration = 10
    sim.compile(verbose=False)
    assert sim._is_compiled
    return {
        "species": set(sim._species_for_sbml.keys()),
        "n_reactions": len(sim._reactions_for_sbml),
    }


def _define_and_compile_model_b() -> dict[str, Any]:
    """Define a different two-species model and compile it."""
    X, Y = BaseSpecies(["X", "Y"])
    _ = X >> Y @ 0.5
    _ = Y >> X @ 0.3
    X(50)
    Y(50)
    sim = Simulation(X | Y)
    sim.duration = 5
    sim.compile(verbose=False)
    assert sim._is_compiled
    return {
        "species": set(sim._species_for_sbml.keys()),
        "n_reactions": len(sim._reactions_for_sbml),
    }


def _define_model_with_inheritance() -> dict[str, Any]:
    """Define a model using inheritance (New) and compile it."""
    Base = BaseSpecies(["Base"])
    Child = New(Base, ["Child"])
    _ = Base >> Zero @ 1
    Base(10)
    Child(20)
    sim = Simulation(Base)
    sim.duration = 5
    sim.compile(verbose=False)
    assert sim._is_compiled
    return {
        "species": set(sim._species_for_sbml.keys()),
        "n_reactions": len(sim._reactions_for_sbml),
    }


def _define_model_with_at_syntax() -> dict[str, Any]:
    """Define a model using @ rate syntax and registry."""
    from mobspy.modules.declarations import get_registry

    A, B = BaseSpecies(["A", "B"])
    A >> B @ 0.5
    B >> A @ 0.3
    A(100)

    reg = get_registry()
    n_registry_reactions = len(reg.reaction_objects)

    sim = Simulation(A | B)
    sim.duration = 5
    sim.compile(verbose=False)
    assert sim._is_compiled
    return {
        "species": set(sim._species_for_sbml.keys()),
        "n_reactions": len(sim._reactions_for_sbml),
        "n_registry": n_registry_reactions,
    }


def _define_model_with_events() -> dict[str, Any]:
    """Define a model using S.at() explicit event API."""
    A, B = BaseSpecies(["A", "B"])
    A >> B @ 0.1
    A(100)

    sim = Simulation(A | B)
    sim.duration = 50
    sim.at(10, {A: 50})
    sim.when(A <= 20, {B: 200})
    sim.compile(verbose=False)
    assert sim._is_compiled
    return {
        "species": set(sim._species_for_sbml.keys()),
        "n_events": len(sim._events_for_sbml),
    }


class TestThreadSafety:
    """Verify that concurrent model definitions don't interfere."""

    def test_two_models_parallel(self) -> None:
        """Two threads define and compile different models simultaneously."""
        results: dict[str, dict[str, Any]] = {}
        errors: list[Exception] = []

        def run_model(name: str, func: Any) -> None:
            try:
                results[name] = func()
            except Exception as e:
                errors.append(e)

        t1 = threading.Thread(target=run_model, args=("a", _define_and_compile_model_a))
        t2 = threading.Thread(target=run_model, args=("b", _define_and_compile_model_b))

        t1.start()
        t2.start()
        t1.join(timeout=30)
        t2.join(timeout=30)

        assert not errors, f"Thread errors: {errors}"
        assert "A" in results["a"]["species"]
        assert results["a"]["n_reactions"] == 2
        assert "X" in results["b"]["species"]
        assert "Y" in results["b"]["species"]
        assert results["b"]["n_reactions"] == 2

    def test_many_threads(self) -> None:
        """Run 8 threads each compiling a model to stress-test."""
        results: list[dict[str, Any] | None] = [None] * 8
        errors: list[Exception] = []

        def run(idx: int) -> None:
            try:
                results[idx] = _define_and_compile_model_a()
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=run, args=(i,)) for i in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=60)

        assert not errors, f"Thread errors: {errors}"
        for r in results:
            assert r is not None
            assert "A" in r["species"]
            assert r["n_reactions"] == 2

    def test_inheritance_parallel(self) -> None:
        """Inheritance (New) works correctly across threads."""
        results: dict[str, dict[str, Any]] = {}
        errors: list[Exception] = []

        def run(name: str, func: Any) -> None:
            try:
                results[name] = func()
            except Exception as e:
                errors.append(e)

        t1 = threading.Thread(
            target=run, args=("inherit", _define_model_with_inheritance)
        )
        t2 = threading.Thread(target=run, args=("simple", _define_and_compile_model_b))
        t1.start()
        t2.start()
        t1.join(timeout=30)
        t2.join(timeout=30)

        assert not errors, f"Thread errors: {errors}"
        assert "Base" in results["inherit"]["species"]
        assert results["inherit"]["n_reactions"] >= 1
        assert "X" in results["simple"]["species"]

    def test_at_syntax_parallel(self) -> None:
        """@ rate syntax with ModelRegistry works across threads."""
        results: dict[str, dict[str, Any]] = {}
        errors: list[Exception] = []

        def run(name: str, func: Any) -> None:
            try:
                results[name] = func()
            except Exception as e:
                errors.append(e)

        t1 = threading.Thread(target=run, args=("at", _define_model_with_at_syntax))
        t2 = threading.Thread(target=run, args=("simple", _define_and_compile_model_b))
        t1.start()
        t2.start()
        t1.join(timeout=30)
        t2.join(timeout=30)

        assert not errors, f"Thread errors: {errors}"
        assert "A" in results["at"]["species"]
        assert results["at"]["n_reactions"] == 2

    def test_events_parallel(self) -> None:
        """S.at() and S.when() work correctly across threads."""
        results: dict[str, dict[str, Any]] = {}
        errors: list[Exception] = []

        def run(name: str, func: Any) -> None:
            try:
                results[name] = func()
            except Exception as e:
                errors.append(e)

        t1 = threading.Thread(target=run, args=("events", _define_model_with_events))
        t2 = threading.Thread(target=run, args=("at", _define_model_with_at_syntax))
        t1.start()
        t2.start()
        t1.join(timeout=30)
        t2.join(timeout=30)

        assert not errors, f"Thread errors: {errors}"
        assert results["events"]["n_events"] == 2
        assert results["at"]["n_reactions"] == 2
