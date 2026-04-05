"""Complex birth/death tests with multiple cell types, inheritance,
characteristics, function-based rates, and units.

Models a tissue ecosystem with stem cells, differentiated cells, and
immune cells. Rates depend on cell identity (is_a) and state
(characteristics) simultaneously.
"""

from __future__ import annotations

from typing import Any

import mobspy
from mobspy import BaseSpecies, Delta, New, Simulation, u
from tests.conftest import CompiledModelAssertions


class TestMultiCellBirthDeath:
    """Birth/death dynamics with inheritance hierarchy and characteristics."""

    def test_three_cell_types_identity_dependent_death(self) -> None:
        """Three cell types inherit from a common ancestor.

        Death rate depends on cell identity via is_a():
        - StemCell: slow turnover (0.001/h)
        - Neuron: moderate (0.01/h)
        - Fibroblast: fast (0.1/h)
        """
        Cell = BaseSpecies()
        StemCell = New(Cell)
        Neuron = New(Cell)
        Fibroblast = New(Cell)

        def death_rate(r: Any) -> Any:
            if r.is_a(StemCell):
                return 0.001 / u.hour * r
            if r.is_a(Neuron):
                return 0.01 / u.hour * r
            return 0.1 / u.hour * r

        Cell >> mobspy.Zero @ death_rate
        StemCell(1000), Neuron(500), Fibroblast(200)
        S = Simulation(StemCell | Neuron | Fibroblast)
        S.level = -1
        S.compile()

        m = CompiledModelAssertions(S)
        m.has_species("StemCell", "Neuron", "Fibroblast")
        m.species_count("StemCell", 1000)
        m.species_count("Neuron", 500)
        m.species_count("Fibroblast", 200)
        m.has_n_reactions(3)
        m.kinetics_contains("StemCell")
        m.kinetics_contains("Neuron")
        m.kinetics_contains("Fibroblast")

    def test_identity_and_characteristic_combined_birth_death(self) -> None:
        """Cells have active/quiescent states. Death rate depends on both
        identity and characteristic.

        - Active StemCells: very low death (0.0005/h)
        - Quiescent StemCells: low death (0.002/h)
        - Active DifferentiatedCells: moderate death (0.05/h)
        - Quiescent DifferentiatedCells: low death (0.01/h)
        """
        Cell = BaseSpecies()
        Cell.active, Cell.quiescent
        StemCell = New(Cell)
        Differentiated = New(Cell)

        def death_rate(r: Any) -> Any:
            if r.is_a(StemCell) and r.active:
                return 0.0005 / u.hour * r
            if r.is_a(StemCell) and r.quiescent:
                return 0.002 / u.hour * r
            if r.active:
                return 0.05 / u.hour * r
            return 0.01 / u.hour * r

        Cell >> mobspy.Zero @ death_rate

        StemCell.active(500), StemCell.quiescent(300)
        Differentiated.active(200), Differentiated.quiescent(100)

        S = Simulation(StemCell | Differentiated)
        S.level = -1
        S.compile()

        m = CompiledModelAssertions(S)
        m.has_species("StemCell.active", "StemCell.quiescent")
        m.has_species("Differentiated.active", "Differentiated.quiescent")
        m.species_count("StemCell.active", 500)
        m.species_count("Differentiated.quiescent", 100)
        # 2 types x 2 states = 4 death reactions
        m.has_n_reactions(4)
        m.kinetics_contains("StemCell.active")
        m.kinetics_contains("Differentiated.quiescent")

    def test_product_type_birth_with_characteristic_gating(self) -> None:
        """Product type Cell = CyclePhase * Phenotype.

        Birth only happens for cells in G2 phase, death for all.
        Rate depends on phenotype via characteristic query.
        """
        CyclePhase, Phenotype = BaseSpecies()
        CyclePhase.G1, CyclePhase.S, CyclePhase.G2
        Phenotype.epithelial, Phenotype.mesenchymal

        Cell = (CyclePhase * Phenotype).named("Cell")

        # Birth: G2 cells divide back into G1. Epithelial divide faster.
        Cell.G2 >> (Cell + Cell.G1) @ (
            lambda r: 0.5 / u.hour if r.epithelial else 0.1 / u.hour
        )

        # Death: all cells, mesenchymal die faster
        Cell >> mobspy.Zero @ (
            lambda r: 0.02 / u.hour if r.mesenchymal else 0.005 / u.hour
        )

        # Phase transitions
        Cell.G1 >> Cell.S @ (0.1 / u.hour)
        Cell.S >> Cell.G2 @ (0.08 / u.hour)

        Cell.G1.epithelial(100), Cell.G1.mesenchymal(50)
        S = Simulation(Cell)
        S.level = -1
        S.compile()

        m = CompiledModelAssertions(S)
        # 3 phases x 2 phenotypes = 6 concrete species
        assert len(m.clean_species) == 6
        m.has_species("Cell.G1.epithelial", "Cell.G2.mesenchymal", "Cell.S.epithelial")
        m.species_count("Cell.G1.epithelial", 100)
        m.species_count("Cell.G1.mesenchymal", 50)
        # Death: 6, birth from G2: 2, phase G1->S: 2, phase S->G2: 2 = 12
        m.has_n_reactions(12)

    def test_second_order_killing_with_identity_check(self) -> None:
        """Immune cells kill target cells. Kill rate depends on both
        the killer identity and the target state.

        - NKCell kills infected targets fast
        - TCell kills infected targets slower
        - Neither kills healthy targets (rate ~ 0)
        """
        Immune = BaseSpecies()
        NKCell = New(Immune)
        TCell = New(Immune)

        Target = BaseSpecies()
        Target.healthy, Target.infected

        def kill_rate(killer: Any, target: Any) -> Any:
            if not target.infected:
                return 1e-15 * u.liter / u.second
            if killer.is_a(NKCell):
                return 1e-7 * u.liter / u.second
            return 1e-8 * u.liter / u.second

        Immune + Target >> Immune @ kill_rate
        NKCell(100), TCell(200)
        Target.healthy(5000), Target.infected(500)

        S = Simulation(NKCell | TCell | Target)
        S.volume = 1 * u.liter
        S.level = -1
        S.compile()

        m = CompiledModelAssertions(S)
        m.has_species("NKCell", "TCell", "Target.healthy", "Target.infected")
        m.species_count("Target.infected", 500)
        # 2 killers x 2 target states = 4 reactions
        m.has_n_reactions(4)
        m.has_reaction_involving("NKCell")
        m.has_reaction_involving("Target.infected")

    def test_three_level_hierarchy_birth_and_death(self) -> None:
        """Three-level hierarchy: Progenitor -> Precursor -> Mature.

        Birth rate decreases down the hierarchy, death increases.
        Uses is_a at each level.
        """
        Cell = BaseSpecies()
        Progenitor = New(Cell)
        Precursor = New(Progenitor)
        Mature = New(Precursor)

        def birth_rate(r: Any) -> Any:
            if r.is_a(Mature):
                return 0.001 / u.hour
            if r.is_a(Precursor):
                return 0.05 / u.hour
            return 0.2 / u.hour

        Cell >> 2 * Cell @ birth_rate

        def death_fn(r: Any) -> Any:
            if r.is_a(Mature):
                return 0.1 / u.hour * r
            if r.is_a(Precursor):
                return 0.02 / u.hour * r
            return 0.005 / u.hour * r

        Cell >> mobspy.Zero @ death_fn

        Progenitor(100), Precursor(500), Mature(2000)
        S = Simulation(Progenitor | Precursor | Mature)
        S.level = -1
        S.compile()

        m = CompiledModelAssertions(S)
        m.has_species("Progenitor", "Precursor", "Mature")
        m.species_count("Progenitor", 100)
        m.species_count("Mature", 2000)
        # 3 birth + 3 death = 6 reactions
        m.has_n_reactions(6)
        m.kinetics_contains("Progenitor")
        m.kinetics_contains("Mature")

    def test_full_tissue_model_with_product_and_hierarchy(self) -> None:
        """Integration test: tissue model combining product types,
        inheritance, characteristics, units, and function rates.

        Cell = Viability * CyclePhase
        Two concrete types: Epithelial(Cell), Stromal(Cell)
        """
        Viability, CyclePhase = BaseSpecies()
        Viability.alive, Viability.apoptotic
        CyclePhase.G1, CyclePhase.S, CyclePhase.G2M

        Cell = (Viability * CyclePhase).named("Cell")
        Epithelial = New(Cell)
        Stromal = New(Cell)

        # Apoptotic clearance (all types)
        Cell.apoptotic >> mobspy.Zero @ (0.5 / u.hour)

        # Alive -> apoptotic transition (type-dependent)
        Cell.alive >> Cell.apoptotic @ (
            lambda r: 0.001 / u.hour if r.is_a(Epithelial) else 0.005 / u.hour
        )

        # Cell cycle progression (alive cells only, characteristic-gated)
        Cell.alive.G1 >> Cell.alive.S @ (
            lambda r: 0.15 / u.hour if r.is_a(Epithelial) else 0.08 / u.hour
        )
        Cell.alive.S >> Cell.alive.G2M @ (0.1 / u.hour)

        # Division: alive G2M cells produce two G1 daughters
        Cell.alive.G2M >> (Cell + Cell.alive.G1) @ (
            lambda r: 0.3 / u.hour if r.is_a(Epithelial) else 0.1 / u.hour
        )

        Epithelial.alive.G1(200), Epithelial.alive.S(50)
        Stromal.alive.G1(300), Stromal.alive.S(80)

        S = Simulation(Epithelial | Stromal)
        S.level = -1
        S.compile()

        m = CompiledModelAssertions(S)
        # Each type has 2 viability x 3 phases = 6 species, x2 types = 12
        assert len(m.clean_species) == 12
        m.has_species("Epithelial.G1.alive", "Epithelial.G2M.apoptotic")
        m.has_species("Stromal.S.alive", "Stromal.G1.apoptotic")
        m.species_count("Epithelial.G1.alive", 200)
        m.species_count("Stromal.G1.alive", 300)

        # Reactions: apoptotic clearance (6), alive->apoptotic (6),
        # G1->S for alive (2), S->G2M for alive (2), division G2M (2) = 18
        m.has_n_reactions(18)


class TestDeltaEvents:
    """Tests for additive (Delta) event assignments."""

    def test_delta_event_generates_additive_sbml(self) -> None:
        """Delta(100) should produce 'species + 100' in SBML event."""
        A = BaseSpecies()
        A >> mobspy.Zero @ 0.1
        A(50)
        S = Simulation(A)
        S.level = -1
        S.at(10, {A: Delta(100)})
        S.compile()

        cm = S._concrete_model
        event = next(iter(cm.events.values()))
        # The assignment value should be "A + 100" (additive)
        assert len(event.assignments) == 1
        variable, value = event.assignments[0]
        assert "A" in variable
        assert "+ 100" in str(value)

    def test_delta_negative(self) -> None:
        """Delta(-50) should produce 'species + -50' in SBML event."""
        A = BaseSpecies()
        A >> mobspy.Zero @ 0.1
        A(200)
        S = Simulation(A)
        S.level = -1
        S.at(5, {A: Delta(-50)})
        S.compile()

        cm = S._concrete_model
        event = next(iter(cm.events.values()))
        variable, value = event.assignments[0]
        assert "A" in variable
        assert "-50" in str(value)

    def test_delta_with_characteristics(self) -> None:
        """Delta works with species that have characteristics."""
        Cell = BaseSpecies()
        Cell.alive, Cell.dead
        Cell.alive >> Cell.dead @ 0.01
        Cell.alive(100), Cell.dead(0)
        S = Simulation(Cell)
        S.level = -1
        S.at(10, {Cell.alive: Delta(50)})
        S.compile()

        cm = S._concrete_model
        event = next(iter(cm.events.values()))
        assert len(event.assignments) == 1
        variable, value = event.assignments[0]
        assert "alive" in variable
        assert "+ 50" in str(value)
