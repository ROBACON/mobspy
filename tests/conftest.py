"""Shared fixtures and utilities for MobsPy tests."""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

import pytest


@pytest.fixture(autouse=True)
def _clear_registry() -> Generator[None, None, None]:
    """Clear the thread-local registry between tests.

    Not required for correctness -- see test_model_isolation.py
    for proof that models in separate functions are independent
    even with a dirty registry. This fixture just prevents
    unbounded memory growth across the test suite.
    """
    from mobspy.modules.session_context import reset_session

    yield
    reset_session()


# ------------------------------------------------------------------
# Semantic test helpers (alternative to golden-file comparison)
# ------------------------------------------------------------------


class CompiledModelAssertions:
    """Helpers for asserting on compiled model structure.

    Use these instead of golden-file comparison when the test
    cares about semantics (species count, reaction structure)
    rather than exact string output.

    Examples:
        >>> from mobspy import BaseSpecies, Simulation
        >>> A, B = BaseSpecies(['A', 'B'])
        >>> _ = A >> B @ 1.0
        >>> A(10)
        >>> S = Simulation(A | B)
        >>> S.compile(verbose=False)  # doctest: +ELLIPSIS
        ...
        >>> a = CompiledModelAssertions(S)
        >>> a.has_species('A', 'B')
        >>> a.species_count('A', 10)
        >>> a.has_n_reactions(1)
    """

    def __init__(self, simulation: Any) -> None:
        self._sim = simulation
        assert simulation._is_compiled, "Simulation must be compiled first"

    @property
    def species(self) -> dict[str, int | float]:
        """Return the species_for_sbml dict."""
        result: dict[str, int | float] = self._sim._species_for_sbml
        return result

    @property
    def reactions(self) -> dict[str, Any]:
        """Return the reactions_for_sbml dict."""
        result: dict[str, Any] = self._sim._reactions_for_sbml
        return result

    def has_species(self, *names: str) -> None:
        """Assert that all named species exist in the compiled model."""
        for name in names:
            assert name in self.species, (
                f"Species '{name}' not found. Available: {sorted(self.species.keys())}"
            )

    def species_count(self, name: str, expected: int | float) -> None:
        """Assert the initial count of a species."""
        self.has_species(name)
        actual = self.species[name]
        assert actual == expected, (
            f"Species '{name}' count: expected {expected}, got {actual}"
        )

    def has_n_reactions(self, n: int, *, exclude_phantom: bool = True) -> None:
        """Assert the number of reactions."""
        if exclude_phantom:
            rxns = {k: v for k, v in self.reactions.items() if "phantom" not in k}
        else:
            rxns = self.reactions
        assert len(rxns) == n, (
            f"Expected {n} reactions, got {len(rxns)}: {list(rxns.keys())}"
        )

    def has_reaction_involving(self, species_name: str) -> None:
        """Assert at least one reaction involves the named species."""
        for rxn in self.reactions.values():
            for _, spe in rxn.reactants:
                if species_name in spe:
                    return
            for _, spe in rxn.products:
                if species_name in spe:
                    return
        msg = f"No reaction involves '{species_name}'"
        raise AssertionError(msg)

    def kinetics_contains(self, substring: str) -> None:
        """Assert at least one reaction's kinetics contains the substring."""
        for rxn in self.reactions.values():
            if substring in rxn.kinetics:
                return
        msg = (
            f"No reaction kinetics contains '{substring}'. "
            f"Kinetics: {[r.kinetics for r in self.reactions.values()]}"
        )
        raise AssertionError(msg)
