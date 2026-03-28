"""Tests for model compilation: basic models, inheritance, species multiplication."""

from __future__ import annotations

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

from .conftest import compare_model, compare_model_ignore_order


class TestBasicModels:
    def test_model_1(self):
        A, B, C = BaseSpecies()
        A + B >> C[1]
        MySim = Simulation(A | B | C)
        MySim.level = -1
        assert compare_model(MySim.compile(), "model_1.txt")

    def test_model_2(self):
        Carnivore, Herbivore = BaseSpecies()
        Cat, Dog = New(Carnivore, 2)
        Carnivore + Herbivore(1 * u.mol) >> Carnivore[1]
        Cat(1 * u.mol), Dog(1 * u.mol)
        MySim = Simulation(Cat | Dog | Herbivore)
        MySim.level = -1
        MySim.volume = 1 * u.meter**2
        assert compare_model(MySim.compile(), "model_2.txt")

    def test_model_3(self):
        MGMT, Blue_Oyster_Cult, The_Smiths = BaseSpecies(3)
        MGMT.eletric_fell, MGMT.little_dark_age, MGMT.kids
        Blue_Oyster_Cult.burning_for_you >> Blue_Oyster_Cult.reaper[1]
        The_Smiths.stop_me >> The_Smiths.charming_man[1]
        Music = MGMT * Blue_Oyster_Cult * The_Smiths
        MySim = Simulation(Music)
        MySim.level = -1
        assert compare_model(MySim.compile(), "model_3.txt")

    def test_model_4(self):
        Bacteria, Virus = BaseSpecies()
        B1, B2 = New(Bacteria, 2)
        V1, V2 = New(Virus, 2)
        Bacteria.not_infected + Virus >> Bacteria.infected[1]
        MySim = Simulation(B1 | B2 | V1 | V2)
        MySim.level = -1
        assert compare_model(MySim.compile(), "model_4.txt")

    def test_model_5(self):
        A = BaseSpecies()
        B, C = New(A, 2)
        A >> 2 * A[1]
        2 * A >> 3 * A[1]
        MySim = Simulation(B | C)
        MySim.level = -1
        assert compare_model(MySim.compile(), "model_5.txt")

    def test_model_6(self):
        A = BaseSpecies()
        B = New(A)
        C = New(A)
        B.b1, B.b2, C.c1, C.c2
        mobspy.Zero >> 2 * A[1]
        MySim = Simulation(B | C)
        MySim.level = -1
        assert compare_model(MySim.compile(), "model_6.txt")

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
                    >> Protein.c(p)
                    + mRNA.c(m)[lambda pro: f"{beta_m}/(1 + ({pro}/{k})^{n})"]
                )
            for m, p in zip(["m1", "m2", "m3"], ["x1", "x2", "x3"], strict=False):
                mRNA.c(m) >> mRNA.c(m) + Protein.c(p)[beta_p]
            Mortal >> mobspy.Zero[lambda r1: gamma_p if r1.is_a(Protein) else gamma_m]
            mobspy.Zero >> Creator[leaky]
            MySim = Simulation(mRNA | Protein)
            MySim.level = -1
            return MySim.compile()

        assert compare_model(oscillator(), "model_7.txt")


class TestInheritanceAndQueries:
    def test_zero_rate_reactions(self):
        A, B = BaseSpecies(2)
        A.a1, A.a2, B.b1, B.b2
        Combination = A * B
        Combination >> mobspy.Zero[lambda r1: 0 if r1.b2 else 1]
        S = Simulation(Combination)
        S.level = -1
        assert compare_model(S.compile(), "model_12.txt")

    def test_double_rate(self):
        A, B = BaseSpecies(2)
        A.a1, A.a2, B.b1, B.b2

        def rate(r1, ball):
            factor1 = 0.5 if r1.a1 else 1
            factor2 = 0.5 if ball.b1 else 1
            return factor1 * factor2

        A + B >> mobspy.Zero[rate]
        S = Simulation(A | B)
        S.level = -1
        assert compare_model(S.compile(), "model_13.txt")

    def test_single_rate(self):
        A, B = BaseSpecies(2)
        A.a1, A.a2, B.b1, B.b2

        def rate(r1):
            factor1 = 0.5 if r1.a1 else 1
            return factor1

        A + B >> mobspy.Zero[rate]
        S = Simulation(A | B)
        S.level = -1
        assert compare_model(S.compile(), "model_14.txt")

    def test_triple_rate(self):
        A, B = BaseSpecies(2)
        A.a1, A.a2, B.b1, B.b2

        def rate(r1, r2, r3):
            factor1 = 0.5 if r1.a1 else 1
            factor2 = 0.5 if r2.b1 else 1
            factor3 = 0.5 if r3.c1 else 1
            return factor1 * factor2 * factor3

        A + B >> mobspy.Zero[rate]
        S = Simulation(A | B)
        S.level = -1
        assert compare_model(S.compile(), "model_13.txt")

    def test_matching_characteristic_rate(self):
        A = BaseSpecies()
        B, C = New(A)
        A.a1, A.a2, A.a3
        B + C >> mobspy.Zero[lambda r1, r2: 100 if A(r1) == A(r2) else 0]
        S = Simulation(B | C)
        S.level = -1
        assert compare_model(S.compile(), "model_43.txt")

    def test_stack_position(self):
        Cell = BaseSpecies()
        Cell >> 2 * Cell[1]
        A, B, C = New(Cell)
        S = Simulation(A | B | C)
        S.level = -1
        compare_model(S.compile(), "model_17.txt")

        def hi():
            Cell = BaseSpecies()
            Cell >> 2 * Cell[1]
            A, B, C = New(Cell)
            S = Simulation(A | B | C)
            S.level = -1
            compare_model(S.compile(), "model_17.txt")

        hi()

        def hi_inside_hi():
            hi()

        hi_inside_hi()


class TestAllOperator:
    def test_all(self):
        A, B = BaseSpecies()
        A.a1, A.a2
        B.b1, B.b2
        C = A * B
        All[C](100)
        C >> All[C][1]
        S = Simulation(C)
        S.level = -1
        assert compare_model(S.compile(), "model_21.txt")

    def test_2_all(self):
        B = BaseSpecies()
        B.b1, B.b2
        C, D = New(B)
        C.c1, C.c2, D.d1, D.d2
        mobspy.Zero >> All[B.b1][1]
        S = Simulation(C | D)
        S.level = -1
        assert compare_model(S.compile(), "model_22.txt")


class TestSetCounts:
    def test_set_counts(self):
        A, C = BaseSpecies()
        A.a1, A.a2
        B = New(A)
        B.b1, B.b2
        model = set_counts({All["B.a1"]: 100, C: 200 * u.mols, "A.a1": 100, A.a2: 50})
        S = Simulation(model)
        S.level = -1
        assert compare_model(S.compile(), "model_23.txt")

    def test_multiple_simulation_counts(self):
        Age, Color, Size = BaseSpecies()
        Age.young, Age.old
        Color.blue, Color.red
        Size.small, Size.big
        Tree = Age * Color * Size
        Tree(100), Tree.red.big(150), All[Tree](10), All[Tree.big](100)
        S = Simulation(Tree)
        S.level = -1
        assert compare_model(S.compile(), "model_28.txt")

        Tree.reset_quantities()
        model = set_counts({All[Tree]: 30, "Tree.blue.old": 100})
        S = Simulation(model)
        S.level = -1
        assert compare_model(S.compile(), "model_29.txt")

    def test_set_counts_parameters(self):
        A = BaseSpecies()
        a = ModelParameters([1, 2])
        A >> 2 * A[a]
        set_counts({"A": a})
        S = Simulation(A)
        S.level = -1
        assert compare_model(S.compile(), "model_32.txt")


class TestDimensions:
    def test_unit_bi_dimension(self):
        A = BaseSpecies()
        A(5 / u.m**2)
        S = Simulation(A)
        S.volume = 2 * u.m**2
        S.level = -1
        assert compare_model(S.compile(), "model_25.txt")

    def test_bi_dimensional_rates(self):
        Ball, Child, Bacteria = BaseSpecies(3)
        Ball(10 / u.meter**2)
        Child(1 / u.meter**2)
        Bacteria(1 * u.mol)
        Bacteria >> mobspy.Zero[1 * u.mol / u.second]
        Ball + Child + Child >> Ball + Child[1e-3 * (u.meter**4) / u.hour]
        Ball + Child >> Ball[1e-3 * (u.meter**2) / u.hour]
        S = Simulation(Ball | Child | Bacteria)
        S.volume = 2 * u.m**2
        S.level = -1
        assert compare_model(S.compile(), "model_26.txt")

    def test_dimension_in_function_only(self):
        A = BaseSpecies()
        A + A >> 3 * A[lambda: 1 * u.milliliter / u.second]
        A(1)
        S = Simulation(A)
        S.level = -1
        assert compare_model(S.compile(), "model_27.txt")

    def test_dimensionless_count(self):
        a = 100 * u.l / u.l
        A = BaseSpecies()
        A >> mobspy.Zero[1]
        A(a)
        S = Simulation(A)
        S.duration = 10
        S.level = -1
        assert compare_model(S.compile(), "model_54.txt")


class TestEmptyArgAndExpressions:
    def test_empty_arguments(self):
        A, B = BaseSpecies()
        A >> mobspy.Zero[lambda: f"{A}*0.01"]
        S = Simulation(A)
        S.duration = 5
        S.level = -1
        assert compare_model(S.compile(), "model_18.txt")

    def test_initial_expression(self):
        A, B, Hey = BaseSpecies()
        D = New(A)
        A >> 2 * A[lambda r: 1 / u.hour * (1 + 10 / r)]
        (
            A + B
            >> mobspy.Zero[
                lambda r1, r2: (
                    (1 * u.millimolar / u.hour)
                    * (1 + 10 * u.millimolar / r1 + 20 * u.millimolar / r2)
                )
            ]
        )
        Hey >> mobspy.Zero[lambda r: 1 / u.hour * (20 * r + 30 * r + 40 * r)]
        D >> 2 * D[lambda r: 20 / u.hour * r]
        S = Simulation(A | B | Hey | D)
        S.level = -1
        assert compare_model(S.compile(), "model_33.txt")

    def test_more_than_used(self):
        A = BaseSpecies()
        mobspy.Zero >> A[lambda r1: 20]
        S = Simulation(A)
        S.level = -1
        assert compare_model(S.compile(), "model_34.txt")

    def test_conversion_outside(self):
        n_0 = 10
        mu_g = 0.2 / u.hour
        Cell, Lysis, AHL, LuxI = BaseSpecies()
        Cell >> 2 * Cell[lambda cell: mu_g * cell * (n_0 - cell)]
        MySim = Simulation(Cell)
        MySim.level = -1
        assert compare_model(MySim.compile(), "model_36.txt")

    def test_first_characteristic_in_reacting_species(self):
        A = BaseSpecies()
        A.something
        B = New(A)
        for a in [1, 2, 3]:
            mobspy.Zero >> B.something.c("at_" + str(a))[1]
        B(1)
        S = Simulation(B)
        S.level = -1
        assert compare_model(S.compile(), "model_37.txt")


class TestReversibleReactions:
    def test_rev(self):
        from mobspy import Rev

        A, B, C = BaseSpecies()
        Rev[A + 4 * B >> C][1, 2]
        Rev[A + 4 * B >> C][lambda r1, r2: (100 - r1) * (100 - r2), lambda r: r**3]
        S = Simulation(A | B | C)
        S.level = -1
        assert compare_model(S.compile(), "model_53.txt")

    def test_new_reversible_reaction_notation(self):
        A = BaseSpecies()
        k1, k2 = ModelParameters(1, 1)
        A >> mobspy.Zero[1, 1]
        A >> mobspy.Zero[1 / (k1 + k2), k1]
        A >> 2 * A[10]
        S = Simulation(A)
        S.level = -1
        assert compare_model(S.compile(), "model_63.txt")


class TestSpeciesNaming:
    def test_silicon_valley(self):
        A = BaseSpecies()
        A.name("\tA")
        A >> mobspy.Zero[1]
        A(200)
        S = Simulation(A)
        S.plot_data = False
        S.duration = 1
        S.level = -1
        assert compare_model(S.compile(), "model_48.txt")

    def test_assignment_similar_species(self):
        A, R, Raa = BaseSpecies()
        A.assign(R * Raa)
        S = Simulation(A | R | Raa)
        S.level = -1
        assert compare_model(S.compile(), "model_55.txt")

    def test_blocked_names(self):
        try:
            _S0 = BaseSpecies()
            _S0 >> mobspy.Zero[1]
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
        S0 >> mobspy.Zero[1]
        S1 >> mobspy.Zero[1]
        S2 >> mobspy.Zero[1]
        S = Simulation(S0 | S1 | S2)
        S.level = -1
        assert compare_model(S.compile(), "model_56.txt")


class TestSBMLGeneration:
    def test_sbml_generation(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1]
        A(100)
        S = Simulation(A)
        S.level = -1
        text = ""
        for sbml in S.generate_sbml():
            text += sbml
        assert compare_model(text, "model_38.txt")

    def test_multi_sim_sbml(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1]
        A(100)
        S1 = Simulation(A)
        S2 = Simulation(A)
        S = S1 + S2
        S.level = -1
        text = ""
        for sbml in S.generate_sbml():
            text += sbml
        assert compare_model(text, "model_39.txt")

    def test_inline_comment(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1]  # Test comment
        assert True


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
                Tree >> Grass[1]
            with Color.blue:
                Tree >> Grass[1]
                Tree(10)
            Tree(9)
            All[Grass](1)
        with mobspy.Any.young.green:
            Tree + Grass >> Tree + Tree[2]

        S1 = Simulation(Tree | Grass)
        S1.level = -1
        assert compare_model(S1.compile(), "model_40.txt")

        Age, Color, Dense = BaseSpecies()
        Age.old, Age.young
        Color.red, Color.green
        Dense.dense, Dense.sparse
        Tree = Age * Color * Dense
        Grass = Age * Color * Dense

        with Age.old, Dense.sparse:
            with mobspy.Any.red:
                Tree >> Grass[1]
            with Color.blue:
                Tree >> Grass[1]
                Tree(10)
            Tree(9)
            All[Grass](1)

        S2 = Simulation(Tree | Grass)
        S2.level = -1
        assert compare_model(S2.compile(), "model_41.txt")

    def test_with_statement_on_any_and_event(self):
        A = BaseSpecies()
        A.a1, A.a2
        S = Simulation(A)
        S.level = -1
        with mobspy.Any.a2, S.event_condition(A <= 0):
            A(100)
        assert compare_model(S.compile(), "model_42.txt")


class TestParameters:
    def test_parameter_operation_in_rate(self):
        A, B = BaseSpecies()
        a = ModelParameters(0.1)
        A >> mobspy.Zero[a]
        B >> mobspy.Zero[2 * a]
        A(100), B(200)
        S1 = Simulation(A | B)
        S1.level = -1
        assert compare_model(S1.compile(), "model_44.txt")

    def test_parameters_as_initial_values(self):
        L, R = BaseSpecies()
        L_0, R_0 = ModelParameters(100, 200)
        L(L_0), R(R_0)
        S = Simulation(L | R)
        S.level = -1
        assert compare_model(S.compile(), "model_65.txt")

    def test_parameters_in_lambda_expression(self):
        L, R = BaseSpecies()
        kf, kr = ModelParameters(1e-3, 1e-3)
        L.sl_0 + R.sr_0 >> L.sl_1 + R.sr_1[kf, lambda r: kr * r]
        S = Simulation(L | R)
        S.level = -1
        assert compare_model(S.compile(), "model_66.txt")

    def test_update_parameter_through_str(self):
        A = BaseSpecies()
        k1 = ModelParameters(0.00000001)
        A >> mobspy.Zero[k1]
        S = Simulation(A)
        S.level = -1
        S.compile()
        S.update_model(["k1", 1])
        assert compare_model(S.generate_sbml()[0], "model_58.txt")

    def test_update_multiple_parameters_in_expression(self):
        A = BaseSpecies()
        k1, k2 = ModelParameters(0.00000001, 10)
        A >> mobspy.Zero[k1 / (10 + k2**4)]
        S = Simulation(A)
        S.level = -1
        S.compile()
        S.update_model([k1, 1], [k2, 1])
        assert compare_model_ignore_order(S.generate_sbml()[0], "model_59.txt")

    def test_update_parameter_with_unit(self):
        A = BaseSpecies()
        k1 = ModelParameters(1 / u.h)
        A >> mobspy.Zero[k1]
        S = Simulation(A)
        S.level = -1
        S.compile()
        S.update_model([k1, 1 / u.s])
        assert compare_model(S.generate_sbml()[0], "model_60.txt")

    def test_species_value_modification(self):
        A = BaseSpecies()
        A.a1, A.a2
        B = New(A)
        B.b1, B.b2
        k1 = ModelParameters(1)
        B >> mobspy.Zero[k1]
        B(100), B.b2(100)
        S = Simulation(B)
        S.level = -1
        S.compile()
        S.update_model([B, 200 / u.l], [B.b2, 300 / u.l])
        assert compare_model_ignore_order(S.generate_sbml()[0], "model_61.txt")

    def test_all_value_modification(self):
        A = BaseSpecies()
        A.a1, A.a2
        B = New(A)
        B.b1, B.b2
        k1 = ModelParameters(1)
        B >> mobspy.Zero[k1]
        B(100), B.b2(100)
        S = Simulation(B)
        S.level = -1
        S.compile()
        S.update_model([All[B], 200 / u.l])
        assert compare_model_ignore_order(S.generate_sbml()[0], "model_62.txt")

    def test_convert_back_parameter(self):
        p = ModelParameters(2 * u.mol / u.l)
        p.convert_to_original_unit()
        assert p.value.magnitude == (2 * u.mol / u.l).magnitude
        assert p.value.units == (2 * u.mol / u.l).units

    def test_update_parameter_for_multi_model(self):
        A, B, C, D = BaseSpecies()
        k1 = ModelParameters([1, 2, 3])
        A >> mobspy.Zero[2 * k1]
        A(100)
        S1 = Simulation(A)
        S1.level = -1
        S1.duration = 10

        A.reset_reactions()
        B >> mobspy.Zero[1]
        B(200)
        S2 = Simulation(A | B)
        S2.level = -1
        S2.duration = 5

        S = S1 + S2
        S.level = -1
        S.compile()
        S.update_model([k1, 1])
        sbml = S.generate_sbml()[0] + "\n" + S.generate_sbml()[1]
        assert compare_model(sbml, "model_57.txt")

    def test_2D_reaction_with_units(self):
        Color, Location = BaseSpecies()
        Color.red, Color.blue
        Location.here, Location.there
        Something = Color * Location
        rate = (  # noqa: E731
            lambda r1, r2: (
                1 * u.decimeter**2 / u.h
                if Location(r1) == Location(r2)
                else 0.5 * u.decimeter**2 / u.h
            )
        )
        2 * Something >> 3 * Something[rate]
        S = Simulation(Something)
        S.level = -1
        S.volume = 1 * u.m**2
        assert compare_model(S.compile(), "model_64.txt")
