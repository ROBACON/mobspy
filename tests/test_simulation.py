"""Tests for simulations (deterministic, stochastic, hybrid, concatenated)."""

from __future__ import annotations

import os

import pytest

import mobspy
from mobspy import (
    All,
    BaseSpecies,
    ModelParameters,
    New,
    Simulation,
    set_counts,
    simlog,
    u,
)
from mobspy.exceptions import MobsPyError

from .conftest import compare_model


@pytest.mark.slow
class TestSimulationExecution:
    def test_average_value(self):
        E = BaseSpecies(1)
        mobspy.Zero >> E[12]
        E >> mobspy.Zero[25]
        MySim = Simulation(E)
        MySim.save_data = False
        MySim.run(plot_data=False)

    def test_hybrid_sim(self):
        A, B = BaseSpecies(2)
        A >> 2 * A[1]
        A(1), B(10)
        S1 = Simulation(A)
        S1.save_data = False
        S1.plot_data = False
        S1.duration = 2

        A.reset_reactions()
        A + B >> mobspy.Zero[0.1]

        S2 = Simulation(A | B)
        S2.method = "stochastic"
        S2.duration = (A <= 0) | (B <= 0)
        S2.level = -1
        S2.plot_data = False

        Sim = S1 + S2
        Sim.run(plot_data=False)

        assert compare_model(Sim.compile(), "model_8.txt")
        assert Sim.fres[A][-1] == 0 or Sim.fres[B][-1] == 0

    def test_concatenated_simulation(self):
        A, B, C = BaseSpecies(3)
        A >> mobspy.Zero[1]
        A(50)
        S1 = Simulation(A)
        S1.plot_data = False
        S1.duration = 5

        B >> mobspy.Zero[1]
        B(50)
        S2 = Simulation(B)
        S2.duration = 5
        S2.plot_data = False

        C >> mobspy.Zero[1]
        C(50)
        S3 = Simulation(C)
        S3.duration = 5
        S3.level = -1

        S = S1 + S2 + S3
        S.run(plot_data=False)
        assert S.fres[A][-1] < 1 and S.fres[B][-1] < 1 and S.fres[C][-1] < 1

    def test_stochastic_event_duration(self):
        A, B = BaseSpecies(2)
        A + B >> mobspy.Zero[0.1]
        A(20), B(20)
        S1 = Simulation(A | B)
        S1.save_data = False
        S1.plot_data = False
        S1.method = "stochastic"
        S1.duration = (A <= 0) | (B <= 0)
        S1.level = -1
        S1.run(plot_data=False)
        R = S1.fres
        assert R[A][0] > 0 and R[B][0] > 0 and R[A][-1] == 0 and R[B][-1] == 0

    def test_reaction_deactivation(self):
        A, R = BaseSpecies(2)
        A + R >> 2 * A + R[1]
        R >> mobspy.Zero[1e-100]
        A(1), R(1)
        S1 = Simulation(A | R)
        S1.level = -1
        S1.duration = 1
        S1.plot_data = False

        S2 = Simulation(A | R)
        S2.duration = 1
        S2.level = -1
        with S2.event_time(0):
            R(0)

        Sim = S1 + S2
        Sim.run(plot_data=False)
        assert (
            Sim.fres[A][0] < Sim.fres[A][-1]
            and Sim.fres[R][0] == 1
            and Sim.fres[R][-1] == 0
        )

    def test_count_assignment(self):
        A = BaseSpecies(1)
        B = New(A)
        A.a1, A.a2
        B.b1 >> mobspy.Zero[1]
        A.a1(100), A.a2(100)
        B.b1(100), B.b2(100)
        S = Simulation(A | B)
        S.level = -1
        S.plot_data = False
        S.duration = 5
        S.run(plot_data=False)
        assert (
            compare_model(S.compile(), "model_11.txt")
            and 150 > S.fres[B][-1] > 100
            and S.fres[A][-1] == 200
        )

    def test_one_value_concatenation_sim(self):
        A, B = BaseSpecies()
        B(200)
        S2 = Simulation(A | B)
        S2.plot_data = False
        S2.level = -1
        S2.duration = 5
        S2.step_size = 1
        S2.duration = (A <= 0) | (B <= 0)
        S2.run(plot_data=False)
        assert len(S2.fres[A]) == 1

    def test_volume_after_sim(self):
        A = BaseSpecies()
        mobspy.Zero >> A[42 * 1 / (u.s * u.milliliter)]
        A >> mobspy.Zero[1]
        S = Simulation(A)
        S.plot_data = False
        S.output_concentration = False
        S.level = -1
        S.volume = 1 * u.milliliter
        S.run(plot_data=False)
        assert int(S.fres[A][-1]) == 42

    def test_changes_after_compilation(self):
        A, B = BaseSpecies()
        A + B >> mobspy.Zero[1]
        A(200), B(200)
        Sim = Simulation(A | B)
        Sim.level = -1
        Sim.compile()
        Sim.duration = 30 * u.hour
        Sim.volume = 1 * u.m**3
        Sim.run(plot_data=False)
        assert Sim._parameters_for_sbml["volume"][0] > 100
        assert Sim.fres["Time"][-1] > 100

    def test_duration_with_run(self):
        A, B = BaseSpecies()
        A + B >> mobspy.Zero[0.01]
        A(10), B(5)
        S = Simulation(A | B)
        S.method = "stochastic"
        S.duration = (A <= 0) | (B <= 0)
        S.run(level=-1, plot_data=False)
        assert S.fres[B][-1] == 0


@pytest.mark.slow
class TestRunArguments:
    def test_run_args(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1]
        A(100)
        S = Simulation(A)
        S.run(duration=1, volume=10, plot_data=False, level=-1, step_size=0.25, jobs=2)
        assert S.__dict__["parameters"]["duration"] == 1
        assert S.__dict__["parameters"]["plot_data"] is False
        assert S.__dict__["parameters"]["step_size"] == 0.25
        assert S.__dict__["parameters"]["level"] == -1
        assert S.__dict__["parameters"]["jobs"] == 2

    def test_unit_args(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1 / u.year]
        A(1 * u.mol)
        S = Simulation(A)
        S.level = -1
        S.run(
            duration=10 * u.year,
            step_size=1 * u.year,
            unit_x=u.year,
            unit_y=u.mol,
            jobs=1,
            level=-1,
            plot_data=False,
        )
        assert S.fres["Time"][-1] > 9.9
        assert S.fres["Time"][1] > 0.99
        assert str(S.__dict__["parameters"]["unit_x"]) == str(1 * u.year)
        assert str(S.__dict__["parameters"]["unit_y"]) == str(1 * u.mol)

    def test_multi_parameters_in_run(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1]
        A(50)
        S1 = Simulation(A)
        S2 = Simulation(A)
        S = S1 + S2
        S.run(
            duration=[2, 3],
            simulation_method=["deterministic", "stochastic"],
            level=-1,
            plot_data=False,
        )
        assert S1.__dict__["parameters"]["simulation_method"] == "deterministic"
        assert S2.__dict__["parameters"]["simulation_method"] == "stochastic"
        assert S1.__dict__["parameters"]["duration"] == 2
        assert S2.__dict__["parameters"]["duration"] == 3

    def test_output_concentration_in_multi_sim(self):
        A, B = BaseSpecies()
        A + B >> mobspy.Zero[0.001]
        A(100), B(200)
        S1 = Simulation(A | B)
        S1.duration = 5 * u.seconds
        S1.volume = 5

        S2 = Simulation(A | B)
        S2.duration = 5
        S2.volume = 100
        S = S1 + S2
        S.level = -1
        S.plot_data = False
        S.output_concentration = True
        S.run(plot_data=False)
        assert S.fres[A][-1] < 10
        assert S.fres[B][-1] < 10

    def test_unit_x_conversion(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1 / u.h]
        A(100)
        S = Simulation(A)
        S.level = -1
        S.step_size = 0.1 * u.h
        S.duration = 1 * u.h
        S.unit_x = u.h
        S.run(plot_data=False)
        assert round(S.fres["Time"][-1]) == 1


class TestMultiParameterSimulation:
    @pytest.mark.slow
    def test_multi_parameter_with_expression(self):
        A = BaseSpecies()
        p = ModelParameters([0.5, 1, 1.5])
        A >> mobspy.Zero[2 * p]
        A(100)
        S = Simulation(A)
        S.run(duration=1, plot_data=False, level=-1)
        assert int(S.results[A][0][-1]) == 36
        assert int(S.results[A][1][-1]) == 13
        assert int(S.results[A][2][-1]) == 4

    @pytest.mark.slow
    def test_double_parameters_with_units(self):
        A = BaseSpecies()
        p1, p2 = ModelParameters([1], [1 / u.hour, 2 / u.hour, 3 / u.hour])
        A >> mobspy.Zero[p1 * p2]
        A(100)
        S = Simulation(A)
        S.run(duration=5 * u.hour, plot_data=False, level=-1)
        assert compare_model(str(S.results), "model_45.txt")

    def test_parameters_with_units(self):
        A = BaseSpecies()
        p = ModelParameters([1 / u.hour, 2 / u.hour, 3 / u.hour])
        A >> mobspy.Zero[p]
        A(100)
        S2 = Simulation(A)
        S2.level = -1
        S2.compile()

    def test_multi_methods_plot(self):
        A = BaseSpecies()
        S1 = Simulation(A)
        S2 = Simulation(A)
        S2.method = "stochastic"
        S = S1 + S2
        S.level = -1
        S.repetitions = 10
        S.compile()
        assert S2.__dict__["parameters"]["plot_type"] == "stochastic"


@pytest.mark.slow
class TestPlotting:
    def test_plotting(self):
        Color, Disease = BaseSpecies()
        Color.blue, Color.red, Color.yellow
        Disease.not_sick, Disease.sick
        Disease.not_sick >> Disease.sick[1]
        Tree = Color * Disease
        Tree.yellow(20), Tree.red(20), Tree.blue(20)
        S = Simulation(Tree)
        S.level = -1
        S.method = "stochastic"
        S.plot_data = False
        S.repetitions = 1
        S.step_size = 0.25
        S.duration = 1
        S.run(plot_data=False)
        S.plot_config.save_to = "tests/plot_output/stochastic_tree.png"
        S.plot_stochastic(Tree.not_sick, Tree.sick)
        S.plot_config.save_to = "tests/plot_output/deterministic_tree.png"
        S.plot(Tree.not_sick, Tree.sick)
        S.plot_config.save_to = "tests/plot_output/constant_tree.png"
        S.plot()
        assert os.path.exists("tests/plot_output/stochastic_tree.png")  # noqa: PTH110
        assert os.path.exists("tests/plot_output/deterministic_tree.png")  # noqa: PTH110
        assert os.path.exists("tests/plot_output/constant_tree.png")  # noqa: PTH110


class TestErrorHandling:
    def test_orthogonal_spaces(self):
        try:
            A, B = BaseSpecies(2)
            A.a, A.b
            C = New(B)
            C.a, C.b
            MySim = Simulation(A | C)
            MySim.level = -1
            MySim.compile()
            assert False
        except (SystemExit, MobsPyError):
            assert True

    def test_dimensional_inconsistency(self):
        try:
            A, B, C = BaseSpecies(3)
            A(1 * u.mol / u.meter**3) + B(1 * u.mol / u.meter**2) >> C[1]
            MySim = Simulation(A | B | C)
            MySim.level = -1
            MySim.compile()
            assert False
        except (SystemExit, MobsPyError):
            assert True

    def test_error_mult(self):
        try:
            D = BaseSpecies(1)
            A, B, C = D * BaseSpecies(3)
            simlog.global_simlog_level = -1
            assert False
        except (SystemExit, MobsPyError):
            assert True

    @pytest.mark.slow
    def test_crash_after_modification(self):
        try:
            A = BaseSpecies()
            S1 = Simulation(A)
            A = BaseSpecies()
            A.a1, A.a2
            S2 = Simulation(A)
            S = S1 + S2
            S.level = -1
            S.run(plot_data=False)
            assert False
        except (SystemExit, MobsPyError):
            assert True

    def test_wrong_dimension_error(self):
        # First case: 1/hour * (1 + 10/dm³/r) is valid in concentration mode
        # because 10/dm³/r simplifies to dimensionless when r is concentration
        A, B = BaseSpecies()
        A >> 2 * A[lambda r: 1 / u.hour * (1 + 10 / u.decimeter**3 / r)]
        S = Simulation(A)
        S.level = -1
        S.compile()

        # Second case: 1/(hour*dm³) * (1 + 10/r) - genuinely wrong dimensions
        try:
            A, B = BaseSpecies()
            A >> 2 * A[lambda r: (1 / (u.hour * u.decimeter**3)) * (1 + 10 / r)]
            S = Simulation(A)
            S.level = -1
            S.compile()
            assert False
        except (SystemExit, MobsPyError):
            assert True

    def test_wrong_rate(self):
        try:
            Ara, aTc = BaseSpecies()
            Ara >> 2 * Ara[aTc]
            S = Simulation(aTc | Ara)
            S.compile()
            assert False
        except (SystemExit, MobsPyError):
            assert True

    @pytest.mark.slow
    def test_shared_parameter_name(self):
        try:
            A = BaseSpecies()
            a = ModelParameters([1, 2])
            a.rename("A")
            A >> 2 * A[a]
            set_counts({"A": a})
            S = Simulation(A)
            S.level = -1
            S.plot_data = False
            S.run(plot_data=False)
            assert False
        except (SystemExit, Exception):
            assert True

    def test_repeated_parameters(self):
        try:
            A = BaseSpecies()
            A.a1, A.a2
            a = ModelParameters([1, 2])
            A >> 2 * A[a]
            All[A](1)
            S1 = Simulation(A)
            S1.duration = 3

            B = BaseSpecies()
            a = ModelParameters([3, 4])
            B >> 2 * B[a]
            B(1)
            S2 = Simulation(A | B)
            S2.duration = 2

            S = S1 + S2
            S.plot_data = False
            S.level = -1
            S.compile()
            assert False
        except (SystemExit, MobsPyError):
            assert True

    def test_proper_unit_context_exit(self):
        _duration = 40
        rate = 1
        init_res = 10000
        init_bact = 1000
        init_atp = 0

        Res, Bact, ATP = BaseSpecies()
        Res(init_res / u.ul)
        Bact(init_bact / u.ul)
        ATP(init_atp / u.ul)
        Res + Bact >> Bact + Bact + ATP[rate * u.ul / u.hours]
        S = Simulation(Res | Bact | ATP)
        S.level = -1
        S.compile()

        try:
            Res, Bact, ATP = BaseSpecies()
            Res(init_res / u.ul)
            Bact(init_bact / u.ul)
            ATP(init_atp / u.ul)
            Res + Bact >> Bact + Bact + ATP[rate * u.ul / u.meters]
            S = Simulation(Res | Bact | ATP)
            S.level = -1
            S.compile()
            assert False
        except (SystemExit, Exception):
            pass

        Res, Bact, ATP = BaseSpecies()
        Res(init_res / u.ul)
        Bact(init_bact / u.ul)
        ATP(init_atp / u.ul)
        Res + Bact >> Bact + Bact + ATP[rate * u.ul / u.hours]
        S = Simulation(Res | Bact | ATP)
        S.level = -1
        S.compile()
        assert True


class TestModelReference:
    def test_model_reference(self):
        Mortal = BaseSpecies()
        A, B = New(Mortal)
        A(100), B(200)
        S1 = Simulation(A | B)
        r1 = [str(spe) for spe in S1.model]
        assert sorted(r1) == ["A", "B"]

        S1.duration = 0.5
        C = New(Mortal)
        C(50)
        S2 = Simulation(A | B | C)
        r1 = [str(spe) for spe in S1.model]
        r2 = [str(spe) for spe in S2.model]
        assert sorted(r1) == ["A", "B"]
        assert sorted(r2) == ["A", "B", "C"]

    @pytest.mark.slow
    def test_replacing_species_name_in_expression(self):
        Resource, R = BaseSpecies()
        death_rate = lambda r1, r2: r1 * r2 * (u.l / u.s)
        Resource + R >> mobspy.Zero[death_rate]
        S = Simulation(Resource | R)
        S.duration = 10
        S.step_size = 5
        S.plot_data = False
        S.level = -1
        S.run(plot_data=False)
        assert True


class TestAntimony:
    def test_antimony_model(self):
        A, TestSpe = BaseSpecies()
        a = ModelParameters([1, 2])
        A + TestSpe >> mobspy.Zero[a]
        A >> 2 * A[0.01]
        A(2), TestSpe(1)
        S = Simulation(A | TestSpe)
        S.level = -1
        S.duration = 10
        S.generate_antimony()[0]

    def test_antimony_compose_model_gen(self):
        A, B, C = BaseSpecies()
        a = ModelParameters(1)
        A >> 2 * A[a]
        A(1)
        S1 = Simulation(A | C)
        S1.duration = 2
        with S1.event_time(1):
            A(10)
        A >> mobspy.Zero[1]
        B >> 2 * B[1e-20]
        B(10)
        S2 = Simulation(A | B | C)
        S2.duration = 5
        S = S1 + S2
        S.level = -1
        S.generate_antimony(compose=True, model_name="test_compose")[0]
