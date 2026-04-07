"""Tests for model compilation: basic models, inheritance, species multiplication."""

from __future__ import annotations

import pytest

import mobspy
from mobspy import (
    All,
    BaseSpecies,
    ModelParameters,
    New,
    Simulation,
    set_counts,
    u,
)
from mobspy.exceptions import MobsPyError


@pytest.mark.compilation
class TestBasicModels:
    def test_model_1(self):
        A, B, C = BaseSpecies()
        A + B >> C @ 1
        MySim = Simulation(A | B | C)
        MySim.level = -1
        MySim.compile(verbose=False)
        cm = MySim._concrete_model

        assert len(cm.species) == 3
        assert cm.species["A"] == 0
        assert cm.species["B"] == 0
        assert cm.species["C"] == 0
        assert "volume" in cm.parameters
        assert cm.parameters["volume"][0] == 1
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 1
        rxn = next(iter(rxns.values()))
        assert rxn.reactants == [(1, "A"), (1, "B")]
        assert rxn.products == [(1, "C")]
        assert "A * B * 1 * c1" in rxn.kinetics

    def test_model_2(self):
        Carnivore, Herbivore = BaseSpecies()
        Cat, Dog = New(Carnivore, 2)
        Carnivore + Herbivore(1 * u.mol) >> Carnivore @ 1
        Cat(1 * u.mol), Dog(1 * u.mol)
        MySim = Simulation(Cat | Dog | Herbivore)
        MySim.level = -1
        MySim.volume = 1 * u.meter**2
        MySim.compile(verbose=False)
        cm = MySim._concrete_model

        assert len(cm.species) == 3
        assert "Cat" in cm.species
        assert "Dog" in cm.species
        assert "Herbivore" in cm.species
        assert cm.species["Cat"] == 1
        assert cm.species["Dog"] == 1
        assert cm.species["Herbivore"] == 1
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 2
        kinetics = [r.kinetics for r in rxns.values()]
        assert any("Cat" in k and "Herbivore" in k for k in kinetics)
        assert any("Dog" in k and "Herbivore" in k for k in kinetics)

    def test_model_3(self):
        MGMT, Blue_Oyster_Cult, The_Smiths = BaseSpecies(3)
        MGMT.eletric_fell, MGMT.little_dark_age, MGMT.kids
        Blue_Oyster_Cult.burning_for_you >> Blue_Oyster_Cult.reaper @ 1
        The_Smiths.stop_me >> The_Smiths.charming_man @ 1
        Music = MGMT * Blue_Oyster_Cult * The_Smiths
        MySim = Simulation(Music)
        MySim.level = -1
        MySim.compile(verbose=False)
        cm = MySim._concrete_model

        assert len(cm.species) == 12
        assert "Music.burning_for_you.charming_man.eletric_fell" in cm.species
        assert "Music.reaper.stop_me.kids" in cm.species
        for v in cm.species.values():
            assert v == 0
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 12

    def test_model_4(self):
        Bacteria, Virus = BaseSpecies()
        B1, B2 = New(Bacteria, 2)
        V1, V2 = New(Virus, 2)
        Bacteria.not_infected + Virus >> Bacteria.infected @ 1
        MySim = Simulation(B1 | B2 | V1 | V2)
        MySim.level = -1
        MySim.compile(verbose=False)
        cm = MySim._concrete_model

        assert len(cm.species) == 6
        assert "B1.infected" in cm.species
        assert "B1.not_infected" in cm.species
        assert "B2.infected" in cm.species
        assert "B2.not_infected" in cm.species
        assert "V1" in cm.species
        assert "V2" in cm.species
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 4
        kinetics = [r.kinetics for r in rxns.values()]
        assert any("B1_dot_not_infected" in k and "V1" in k for k in kinetics)
        assert any("B1_dot_not_infected" in k and "V2" in k for k in kinetics)
        assert any("B2_dot_not_infected" in k and "V1" in k for k in kinetics)
        assert any("B2_dot_not_infected" in k and "V2" in k for k in kinetics)

    def test_model_5(self):
        A = BaseSpecies()
        B, C = New(A, 2)
        A >> 2 * A @ 1
        2 * A >> 3 * A @ 1
        MySim = Simulation(B | C)
        MySim.level = -1
        MySim.compile(verbose=False)
        cm = MySim._concrete_model

        assert len(cm.species) == 2
        assert "B" in cm.species
        assert "C" in cm.species
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 6
        reactant_lists = [r.reactants for r in rxns.values()]
        assert any(r == [(1, "B")] for r in reactant_lists)
        assert any(r == [(1, "C")] for r in reactant_lists)
        assert any(r == [(2, "B")] for r in reactant_lists)
        assert any(r == [(2, "C")] for r in reactant_lists)

    def test_model_6(self):
        A = BaseSpecies()
        B = New(A)
        C = New(A)
        B.b1, B.b2, C.c1, C.c2
        mobspy.Zero >> 2 * A @ 1
        MySim = Simulation(B | C)
        MySim.level = -1
        MySim.compile(verbose=False)
        cm = MySim._concrete_model

        assert len(cm.species) == 4
        assert "B.b1" in cm.species
        assert "B.b2" in cm.species
        assert "C.c1" in cm.species
        assert "C.c2" in cm.species
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 4
        for rxn in rxns.values():
            assert rxn.reactants == []
            assert "c1" in rxn.kinetics

    def test_model_7(self):
        def oscillator(
            beta_m=5, beta_p=10, gamma_m=1, gamma_p=0.01, k=1, n=4, leaky=0.0001
        ):
            Mortal, Creator = BaseSpecies(2)
            mRNA = Mortal * Creator
            Protein = New(Mortal)
            for m, p in zip(["m1", "m2", "m3"], ["x2", "x3", "x1"], strict=False):
                (
                    Protein.c(p)
                    >> (Protein.c(p) + mRNA.c(m))
                    @ (lambda pro: f"{beta_m}/(1 + ({pro}/{k})^{n})")
                )
            for m, p in zip(["m1", "m2", "m3"], ["x1", "x2", "x3"], strict=False):
                mRNA.c(m) >> (mRNA.c(m) + Protein.c(p)) @ beta_p
            Mortal >> mobspy.Zero @ (
                lambda r1: gamma_p if r1.is_a(Protein) else gamma_m
            )
            mobspy.Zero >> Creator @ leaky
            MySim = Simulation(mRNA | Protein)
            MySim.level = -1
            MySim.compile(verbose=False)
            return MySim._concrete_model

        cm = oscillator()
        assert len(cm.species) == 6
        assert "Protein.x1" in cm.species
        assert "Protein.x2" in cm.species
        assert "Protein.x3" in cm.species
        assert "mRNA.m1" in cm.species
        assert "mRNA.m2" in cm.species
        assert "mRNA.m3" in cm.species
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 13
        kinetics = [r.kinetics for r in rxns.values()]
        assert any("5/(1 + (Protein_dot_x" in k for k in kinetics)
        assert any("0.0001 * c1" in k for k in kinetics)
        assert any("* 10" in k for k in kinetics)
        assert any("* 0.01" in k for k in kinetics)
        assert any("* 1" in k for k in kinetics)


@pytest.mark.compilation
class TestInheritanceAndQueries:
    def test_zero_rate_reactions(self):
        A, B = BaseSpecies(2)
        A.a1, A.a2, B.b1, B.b2
        Combination = A * B
        Combination >> mobspy.Zero @ (lambda r1: 0 if r1.b2 else 1)
        S = Simulation(Combination)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 4
        assert "Combination.a1.b1" in cm.species
        assert "Combination.a1.b2" in cm.species
        assert "Combination.a2.b1" in cm.species
        assert "Combination.a2.b2" in cm.species
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 2
        kinetics = sorted(r.kinetics for r in rxns.values())
        assert any("Combination_dot_a1_dot_b1" in k for k in kinetics)
        assert any("Combination_dot_a2_dot_b1" in k for k in kinetics)

    def test_double_rate(self):
        A, B = BaseSpecies(2)
        A.a1, A.a2, B.b1, B.b2

        def rate(r1, ball):
            factor1 = 0.5 if r1.a1 else 1
            factor2 = 0.5 if ball.b1 else 1
            return factor1 * factor2

        A + B >> mobspy.Zero @ rate
        S = Simulation(A | B)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 4
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 4
        kinetics = {r.kinetics for r in rxns.values()}
        assert any("0.25" in k for k in kinetics)
        assert any("0.5" in k for k in kinetics)
        assert any(" 1 * c1" in k for k in kinetics)

    def test_single_rate(self):
        A, B = BaseSpecies(2)
        A.a1, A.a2, B.b1, B.b2

        def rate(r1):
            return 0.5 if r1.a1 else 1

        A + B >> mobspy.Zero @ rate
        S = Simulation(A | B)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 4
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 4
        kinetics = sorted(r.kinetics for r in rxns.values())
        assert any("0.5" in k for k in kinetics)
        assert any(" 1 * c1" in k for k in kinetics)

    def test_triple_rate(self):
        A, B = BaseSpecies(2)
        A.a1, A.a2, B.b1, B.b2

        def rate(r1, r2, r3):
            factor1 = 0.5 if r1.a1 else 1
            factor2 = 0.5 if r2.b1 else 1
            factor3 = 0.5 if r3.c1 else 1
            return factor1 * factor2 * factor3

        A + B >> mobspy.Zero @ rate
        S = Simulation(A | B)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 4
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 4
        kinetics = {r.kinetics for r in rxns.values()}
        assert any("0.25" in k for k in kinetics)
        assert any("0.5" in k for k in kinetics)

    def test_matching_characteristic_rate(self):
        A = BaseSpecies()
        B, C = New(A)
        A.a1, A.a2, A.a3
        B + C >> mobspy.Zero @ (lambda r1, r2: 100 if A(r1) == A(r2) else 0)
        S = Simulation(B | C)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 6
        for name in [
            "B.a1",
            "B.a2",
            "B.a3",
            "C.a1",
            "C.a2",
            "C.a3",
        ]:
            assert name in cm.species
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 3
        kinetics = [r.kinetics for r in rxns.values()]
        assert all("100" in k for k in kinetics)

    def test_stack_position(self):
        Cell = BaseSpecies()
        Cell >> 2 * Cell @ 1
        A, B, C = New(Cell)
        S = Simulation(A | B | C)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 3
        assert "A" in cm.species
        assert "B" in cm.species
        assert "C" in cm.species
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 3
        for rxn in rxns.values():
            assert rxn.products[0][0] == 2

        def hi():
            Cell = BaseSpecies()
            Cell >> 2 * Cell @ 1
            A, B, C = New(Cell)
            S = Simulation(A | B | C)
            S.level = -1
            S.compile(verbose=False)
            cm = S._concrete_model
            assert len(cm.species) == 3
            rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
            assert len(rxns) == 3

        hi()

        def hi_inside_hi():
            hi()

        hi_inside_hi()


@pytest.mark.compilation
class TestAllOperator:
    def test_all(self):
        A, B = BaseSpecies()
        A.a1, A.a2
        B.b1, B.b2
        C = A * B
        All[C](100)
        C >> All[C] @ 1
        S = Simulation(C)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 4
        for name in [
            "C.a1.b1",
            "C.a1.b2",
            "C.a2.b1",
            "C.a2.b2",
        ]:
            assert name in cm.species
            assert cm.species[name] == 100
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 16

    def test_2_all(self):
        B = BaseSpecies()
        B.b1, B.b2
        C, D = New(B)
        C.c1, C.c2, D.d1, D.d2
        mobspy.Zero >> All[B.b1] @ 1
        S = Simulation(C | D)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 8
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 4
        for rxn in rxns.values():
            assert rxn.reactants == []
            assert "b1" in rxn.products[0][1]


@pytest.mark.compilation
class TestSetCounts:
    def test_set_counts(self):
        A, C = BaseSpecies()
        A.a1, A.a2
        B = New(A)
        B.b1, B.b2
        model = set_counts({All["B.a1"]: 100, C: 200 * u.mols, "A.a1": 100, A.a2: 50})
        S = Simulation(model)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert cm.species["A.a1"] == 100
        assert cm.species["A.a2"] == 50
        assert cm.species["B.a1.b1"] == 100
        assert cm.species["B.a1.b2"] == 100
        assert cm.species["C"] == 200
        assert len(cm.species) == 7

    def test_multiple_simulation_counts(self):
        Age, Color, Size = BaseSpecies()
        Age.young, Age.old
        Color.blue, Color.red
        Size.small, Size.big
        Tree = Age * Color * Size
        Tree(100), Tree.red.big(150), All[Tree](10), All[Tree.big](100)
        S = Simulation(Tree)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 8
        assert cm.species["Tree.big.red.young"] == 150
        assert cm.species["Tree.big.blue.old"] == 100
        assert cm.species["Tree.big.blue.young"] == 100
        assert cm.species["Tree.big.red.old"] == 100
        assert cm.species["Tree.small.blue.young"] == 100
        assert cm.species["Tree.small.red.young"] == 10
        assert cm.species["Tree.small.red.old"] == 10
        assert cm.species["Tree.small.blue.old"] == 10

        Tree.reset_quantities()
        model = set_counts({All[Tree]: 30, "Tree.blue.old": 100})
        S = Simulation(model)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 8
        assert cm.species["Tree.small.blue.old"] == 100
        for name, val in cm.species.items():
            if name != "Tree_dot_small_dot_blue_dot_old":
                assert val == 30

    def test_set_counts_parameters(self):
        A = BaseSpecies()
        a = ModelParameters([1, 2])
        A >> 2 * A @ a
        set_counts({"A": a})
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 1
        assert str(cm.species["A"]) == "a" or cm.species["A"] == "a"
        assert "a" in cm.parameters
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 1
        rxn = next(iter(rxns.values()))
        assert "a" in rxn.kinetics


@pytest.mark.compilation
class TestDimensions:
    def test_unit_bi_dimension(self):
        A = BaseSpecies()
        A(5 / u.m**2)
        S = Simulation(A)
        S.volume = 2 * u.m**2
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 1
        assert cm.species["A"] == pytest.approx(10.0)
        assert cm.parameters["volume"][0] == pytest.approx(2.0, rel=1e-6)

    def test_bi_dimensional_rates(self):
        Ball, Child, Bacteria = BaseSpecies(3)
        Ball(10 / u.meter**2)
        Child(1 / u.meter**2)
        Bacteria(1 * u.mol)
        Bacteria >> mobspy.Zero @ (1 * u.mol / u.second)
        Ball + Child + Child >> (Ball + Child) @ (1e-3 * (u.meter**4) / u.hour)
        Ball + Child >> Ball @ (1e-3 * (u.meter**2) / u.hour)
        S = Simulation(Ball | Child | Bacteria)
        S.volume = 2 * u.m**2
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 3
        assert cm.species["Ball"] == pytest.approx(20.0)
        assert cm.species["Child"] == pytest.approx(2.0)
        assert cm.species["Bacteria"] == 1
        assert cm.parameters["volume"][0] == 2
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 3

    def test_dimension_in_function_only(self):
        A = BaseSpecies()
        A + A >> 3 * A @ (lambda: 1 * u.milliliter / u.second)
        A(1)
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 1
        assert cm.species["A"] == 1
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 1
        rxn = next(iter(rxns.values()))
        assert rxn.reactants == [(2, "A")]
        assert rxn.products == [(3, "A")]
        assert "0.001" in rxn.kinetics

    def test_dimensionless_count(self):
        a = 100 * u.l / u.l
        A = BaseSpecies()
        A >> mobspy.Zero @ 1
        A(a)
        S = Simulation(A)
        S.duration = 10
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 1
        assert cm.species["A"] == pytest.approx(100.0)
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 1


@pytest.mark.compilation
class TestEmptyArgAndExpressions:
    def test_empty_arguments(self):
        A, B = BaseSpecies()
        A >> mobspy.Zero @ (lambda: f"{A}*0.01")
        S = Simulation(A)
        S.duration = 5
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 1
        assert cm.species["A"] == 0
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 1
        rxn = next(iter(rxns.values()))
        assert "A" in rxn.kinetics
        assert "0.01" in rxn.kinetics

    def test_initial_expression(self):
        A, B, Hey = BaseSpecies()
        D = New(A)
        A >> 2 * A @ (lambda r: 1 / u.hour * (1 + 10 / r))
        (
            A + B
            >> mobspy.Zero
            @ (
                lambda r1, r2: (
                    (1 * u.millimolar / u.hour)
                    * (1 + 10 * u.millimolar / r1 + 20 * u.millimolar / r2)
                )
            )
        )
        Hey >> mobspy.Zero @ (lambda r: 1 / u.hour * (20 * r + 30 * r + 40 * r))
        D >> 2 * D @ (lambda r: 20 / u.hour * r)
        S = Simulation(A | B | Hey | D)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 4
        for name in ["A", "B", "D", "Hey"]:
            assert name in cm.species
            assert cm.species[name] == 0
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 6

    def test_more_than_used(self):
        A = BaseSpecies()
        mobspy.Zero >> A @ (lambda r1: 20)
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 1
        assert cm.species["A"] == 0
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 1
        rxn = next(iter(rxns.values()))
        assert rxn.reactants == []
        assert rxn.products == [(1, "A")]
        assert "20 * c1" in rxn.kinetics

    def test_conversion_outside(self):
        n_0 = 10
        mu_g = 0.2 / u.hour
        Cell, Lysis, AHL, LuxI = BaseSpecies()
        Cell >> 2 * Cell @ (lambda cell: mu_g * cell * (n_0 - cell))
        MySim = Simulation(Cell)
        MySim.level = -1
        MySim.compile(verbose=False)
        cm = MySim._concrete_model

        assert len(cm.species) == 1
        assert "Cell" in cm.species
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 1
        rxn = next(iter(rxns.values()))
        assert rxn.reactants == [(1, "Cell")]
        assert rxn.products == [(2, "Cell")]
        assert "Cell" in rxn.kinetics
        assert "10" in rxn.kinetics

    def test_first_characteristic_in_reacting_species(self):
        A = BaseSpecies()
        A.something
        B = New(A)
        for a in [1, 2, 3]:
            mobspy.Zero >> B.something.c("at_" + str(a)) @ 1
        B(1)
        S = Simulation(B)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 3
        assert "B.at_1.something" in cm.species
        assert "B.at_2.something" in cm.species
        assert "B.at_3.something" in cm.species
        assert cm.species["B.at_1.something"] == 1
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 3
        for rxn in rxns.values():
            assert rxn.reactants == []
            assert "c1" in rxn.kinetics


@pytest.mark.compilation
class TestReversibleReactions:
    def test_rev(self):
        A, B, C = BaseSpecies()
        A + 4 * B >> C @ (1, 2)
        A + 4 * B >> C @ (lambda r1, r2: (100 - r1) * (100 - r2), lambda r: r**3)
        S = Simulation(A | B | C)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 3
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 4
        kinetics = [r.kinetics for r in rxns.values()]
        assert any("c1" in k for k in kinetics)
        assert any("(100-A)" in k and "(100-B)" in k for k in kinetics)
        assert any("C^3" in k for k in kinetics)
        assert any("C * 2" in k for k in kinetics)

    def test_new_reversible_reaction_notation(self):
        A = BaseSpecies()
        k1, k2 = ModelParameters(1, 1)
        A >> mobspy.Zero @ (1, 1)
        A >> mobspy.Zero @ (1 / (k1 + k2), k1)
        A >> 2 * A @ 10
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 1
        assert cm.species["A"] == 0
        assert "k1" in cm.parameters
        assert "k2" in cm.parameters
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 5
        kinetics = [r.kinetics for r in rxns.values()]
        assert any("10" in k for k in kinetics)
        assert any("1/(k1+k2)" in k for k in kinetics)
        assert any("k1 * c1" in k for k in kinetics)
        reverse_rxns = [r for r in rxns.values() if r.reactants == []]
        assert len(reverse_rxns) == 2


@pytest.mark.compilation
class TestSpeciesNaming:
    def test_silicon_valley(self):
        A = BaseSpecies()
        A.name("\tA")
        A >> mobspy.Zero @ 1
        A(200)
        S = Simulation(A)
        S.plot_data = False
        S.duration = 1
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 1
        assert "A" in cm.species
        assert cm.species["A"] == 200
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 1
        rxn = next(iter(rxns.values()))
        assert "A * 1" in rxn.kinetics

    def test_assignment_similar_species(self):
        A, R, Raa = BaseSpecies()
        A.assign(R * Raa)
        S = Simulation(A | R | Raa)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 3
        assert "A" in cm.species
        assert "R" in cm.species
        assert "Raa" in cm.species
        assert len(cm.assignments) == 1
        asgn = next(iter(cm.assignments.values()))
        assert asgn.species == "A"
        assert "R" in asgn.expression
        assert "Raa" in asgn.expression

    def test_blocked_names(self):
        try:
            _S0 = BaseSpecies()
            _S0 >> mobspy.Zero @ 1
            assert False
        except (SystemExit, MobsPyError):
            assert True

    def test_blocked_names_2(self):
        try:
            _S1 = BaseSpecies()
            assert False
        except (SystemExit, MobsPyError):
            pass

        S0, S1, S2 = BaseSpecies()
        S0 >> mobspy.Zero @ 1
        S1 >> mobspy.Zero @ 1
        S2 >> mobspy.Zero @ 1
        S = Simulation(S0 | S1 | S2)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 3
        assert "S0" in cm.species
        assert "S1" in cm.species
        assert "S2" in cm.species
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 3


@pytest.mark.compilation
class TestSBMLGeneration:
    def test_sbml_generation(self):
        A = BaseSpecies()
        A >> mobspy.Zero @ 1
        A(100)
        S = Simulation(A)
        S.level = -1
        text = ""
        for sbml in S.generate_sbml():
            text += sbml
        assert "<?xml" in text
        assert "<sbml" in text
        assert '<species id="A"' in text
        assert 'initialAmount="100"' in text
        assert "<listOfReactions>" in text
        assert 'id="reaction_0"' in text

    def test_multi_sim_sbml(self):
        A = BaseSpecies()
        A >> mobspy.Zero @ 1
        A(100)
        S1 = Simulation(A)
        S2 = Simulation(A)
        S = S1 + S2
        S.level = -1
        text = ""
        for sbml in S.generate_sbml():
            text += sbml
        assert text.count("<?xml") == 2
        assert text.count("<sbml") == 2
        assert text.count('<species id="A"') == 2

    def test_inline_comment(self):
        A = BaseSpecies()
        A >> mobspy.Zero @ 1  # Test comment
        assert True


@pytest.mark.compilation
class TestWithStatement:
    def test_with_statement_any_and_species_characteristics(self):
        Age, Color, Dense = BaseSpecies()
        Age.old, Age.young
        Color.red, Color.green
        Dense.dense, Dense.sparse
        Tree = Age * Color * Dense
        Grass = Age * Color * Dense

        with Age.old, Dense.sparse:
            with Color.red:
                Tree >> Grass @ 1
            with Color.blue:
                Tree >> Grass @ 1
                Tree(10)
            Tree(9)
            All[Grass](1)
        with mobspy.Any.young.green:
            Tree + Grass >> (Tree + Tree) @ 2

        S1 = Simulation(Tree | Grass)
        S1.level = -1
        S1.compile(verbose=False)
        cm1 = S1._concrete_model

        assert len(cm1.species) == 24
        assert cm1.species["Tree.red.sparse.old"] == 9
        assert cm1.species["Tree.blue.sparse.old"] == 10
        assert cm1.species["Grass.blue.sparse.old"] == 1
        assert cm1.species["Grass.red.sparse.old"] == 1
        assert cm1.species["Grass.green.sparse.old"] == 1
        rxns1 = {k: v for k, v in cm1.reactions.items() if "phantom" not in k}
        assert len(rxns1) == 6
        kinetics1 = [r.kinetics for r in rxns1.values()]
        assert any("Tree_dot_blue_dot_sparse_dot_old" in k for k in kinetics1)
        assert any("Tree_dot_red_dot_sparse_dot_old" in k for k in kinetics1)

        Age, Color, Dense = BaseSpecies()
        Age.old, Age.young
        Color.red, Color.green
        Dense.dense, Dense.sparse
        Tree = Age * Color * Dense
        Grass = Age * Color * Dense

        with Age.old, Dense.sparse:
            with mobspy.Any.red:
                Tree >> Grass @ 1
            with Color.blue:
                Tree >> Grass @ 1
                Tree(10)
            Tree(9)
            All[Grass](1)

        S2 = Simulation(Tree | Grass)
        S2.level = -1
        S2.compile(verbose=False)
        cm2 = S2._concrete_model

        assert len(cm2.species) == 24
        rxns2 = {k: v for k, v in cm2.reactions.items() if "phantom" not in k}
        assert len(rxns2) == 2
        kinetics2 = [r.kinetics for r in rxns2.values()]
        assert any("Tree_dot_blue_dot_sparse_dot_old" in k for k in kinetics2)
        assert any("Tree_dot_red_dot_sparse_dot_old" in k for k in kinetics2)

    def test_with_statement_on_any_and_event(self):
        A = BaseSpecies()
        A.a1, A.a2
        S = Simulation(A)
        S.level = -1
        with mobspy.Any.a2, S.event_condition(A <= 0):
            A(100)
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 2
        assert "A.a1" in cm.species
        assert "A.a2" in cm.species
        assert len(cm.events) == 1
        event = next(iter(cm.events.values()))
        assert "<=" in event.trigger or "le" in event.trigger.lower()
        assert any("A_dot_a2" in str(a) for a in event.assignments)


@pytest.mark.compilation
class TestParameters:
    def test_parameter_operation_in_rate(self):
        A, B = BaseSpecies()
        a = ModelParameters(0.1)
        A >> mobspy.Zero @ a
        B >> mobspy.Zero @ (2 * a)
        A(100), B(200)
        S1 = Simulation(A | B)
        S1.level = -1
        S1.compile(verbose=False)
        cm = S1._concrete_model

        assert cm.species["A"] == 100
        assert cm.species["B"] == 200
        assert "a" in cm.parameters
        assert cm.parameters["a"][0] == pytest.approx(0.1)
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 2
        kinetics = [r.kinetics for r in rxns.values()]
        assert any("A * a" in k for k in kinetics)
        assert any("(2*a)" in k for k in kinetics)

    def test_parameters_as_initial_values(self):
        L, R = BaseSpecies()
        L_0, R_0 = ModelParameters(100, 200)
        L(L_0), R(R_0)
        S = Simulation(L | R)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert str(cm.species["L"]) == "L_0" or cm.species["L"] == "L_0"
        assert str(cm.species["R"]) == "R_0" or cm.species["R"] == "R_0"
        assert "L_0" in cm.parameters
        assert "R_0" in cm.parameters
        assert cm.parameters["L_0"][0] == 100
        assert cm.parameters["R_0"][0] == 200

    def test_parameters_in_lambda_expression(self):
        L, R = BaseSpecies()
        kf, kr = ModelParameters(1e-3, 1e-3)
        L.sl_0 + R.sr_0 >> (L.sl_1 + R.sr_1) @ (kf, lambda r: kr * r)
        S = Simulation(L | R)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 4
        for name in ["L.sl_0", "L.sl_1", "R.sr_0", "R.sr_1"]:
            assert name in cm.species
        assert "kf" in cm.parameters
        assert "kr" in cm.parameters
        assert cm.parameters["kf"][0] == pytest.approx(0.001)
        assert cm.parameters["kr"][0] == pytest.approx(0.001)
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 2
        kinetics = [r.kinetics for r in rxns.values()]
        assert any("kf" in k for k in kinetics)
        assert any("kr" in k for k in kinetics)

    def test_update_parameter_through_str(self):
        A = BaseSpecies()
        k1 = ModelParameters(0.00000001)
        A >> mobspy.Zero @ k1
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)
        S.update_model(["k1", 1])
        sbml = S.generate_sbml()[0]
        assert '<parameter id="k1" value="1"' in sbml
        assert '<species id="A"' in sbml
        assert "<listOfReactions>" in sbml

    def test_update_multiple_parameters_in_expression(self):
        A = BaseSpecies()
        k1, k2 = ModelParameters(0.00000001, 10)
        A >> mobspy.Zero @ (k1 / (10 + k2**4))
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)
        S.update_model([k1, 1], [k2, 1])
        sbml = S.generate_sbml()[0]
        assert '<parameter id="k1" value="1"' in sbml
        assert '<parameter id="k2" value="1"' in sbml
        assert '<species id="A"' in sbml

    def test_update_parameter_with_unit(self):
        A = BaseSpecies()
        k1 = ModelParameters(1 / u.h)
        A >> mobspy.Zero @ k1
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)
        S.update_model([k1, 1 / u.s])
        sbml = S.generate_sbml()[0]
        assert '<parameter id="k1" value="1"' in sbml
        assert '<species id="A"' in sbml

    def test_species_value_modification(self):
        A = BaseSpecies()
        A.a1, A.a2
        B = New(A)
        B.b1, B.b2
        k1 = ModelParameters(1)
        B >> mobspy.Zero @ k1
        B(100), B.b2(100)
        S = Simulation(B)
        S.level = -1
        S.compile(verbose=False)
        S.update_model([B, 200 / u.l], [B.b2, 300 / u.l])
        sbml = S.generate_sbml()[0]
        assert "B_dot_a1_dot_b1" in sbml
        assert "B_dot_a1_dot_b2" in sbml
        assert 'initialAmount="300"' in sbml
        assert 'initialAmount="200"' in sbml

    def test_all_value_modification(self):
        A = BaseSpecies()
        A.a1, A.a2
        B = New(A)
        B.b1, B.b2
        k1 = ModelParameters(1)
        B >> mobspy.Zero @ k1
        B(100), B.b2(100)
        S = Simulation(B)
        S.level = -1
        S.compile(verbose=False)
        S.update_model([All[B], 200 / u.l])
        sbml = S.generate_sbml()[0]
        assert sbml.count('initialAmount="200"') == 4

    def test_convert_back_parameter(self):
        p = ModelParameters(2 * u.mol / u.l)
        p.convert_to_original_unit()
        assert p.value.magnitude == (2 * u.mol / u.l).magnitude
        assert p.value.units == (2 * u.mol / u.l).units

    def test_update_parameter_for_multi_model(self):
        A, B, C, D = BaseSpecies()
        k1 = ModelParameters([1, 2, 3])
        A >> mobspy.Zero @ (2 * k1)
        A(100)
        S1 = Simulation(A)
        S1.level = -1
        S1.duration = 10

        A.reset_reactions()
        B >> mobspy.Zero @ 1
        B(200)
        S2 = Simulation(A | B)
        S2.level = -1
        S2.duration = 5

        S = S1 + S2
        S.level = -1
        S.compile(verbose=False)
        S.update_model([k1, 1])
        sbml0 = S.generate_sbml()[0]
        sbml1 = S.generate_sbml()[1]
        assert '<parameter id="k1" value="1"' in sbml0
        assert '<species id="A"' in sbml0
        assert 'initialAmount="100"' in sbml0
        assert '<species id="A"' in sbml1
        assert '<species id="B"' in sbml1
        assert 'initialAmount="200"' in sbml1

    def test_2D_reaction_with_units(self):
        Color, Location = BaseSpecies()
        Color.red, Color.blue
        Location.here, Location.there
        Something = Color * Location
        rate = lambda r1, r2: (
            1 * u.decimeter**2 / u.h
            if Location(r1) == Location(r2)
            else 0.5 * u.decimeter**2 / u.h
        )
        2 * Something >> 3 * Something @ rate
        S = Simulation(Something)
        S.level = -1
        S.volume = 1 * u.m**2
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 4
        for name in [
            "Something.blue.here",
            "Something.blue.there",
            "Something.red.here",
            "Something.red.there",
        ]:
            assert name in cm.species
        assert cm.parameters["volume"][0] == pytest.approx(1.0, rel=1e-6)
        rxns = {k: v for k, v in cm.reactions.items() if "phantom" not in k}
        assert len(rxns) == 16
