"""Tests for events (time-based, condition-based) and logic operators."""

from __future__ import annotations

import mobspy
from mobspy import All, BaseSpecies, New, Simulation, set_counts, simlog, u
from mobspy.exceptions import MobsPyError

from .conftest import compare_model


class TestTimeEvents:
    def test_event_type(self):
        A, B, C, D, E, F = BaseSpecies(6)
        A + B >> mobspy.Zero[1]
        A(50), B(50), C(0)
        S = Simulation(A | B | C | D | E | F)
        S.plot_data = False
        S.level = -1
        with S.event_time(0):
            F(1)
        with S.event_condition((A <= 1) & (B <= 1)):
            C(1)
        with S.event_condition((A <= 1) & (B <= 1)):
            D(1)
        with S.event_condition(B <= 1):
            E(1)
        S.duration = 5
        assert compare_model(S.compile(), "model_15.txt")

    def test_reacting_species_event(self):
        B = BaseSpecies(1)
        B.b1, B.b2
        A = New(B)
        B >> mobspy.Zero[1]
        A.b2 >> mobspy.Zero[0.5]
        A.a1 >> mobspy.Zero[1]
        A.a1(100), A.b2(100), B.b1(100)
        S = Simulation(A | B)
        S.level = -1
        with S.event_condition((A.a1 <= 10) & (B.b1 <= 10)):
            A.a1(100)
        S.duration = 5
        assert compare_model(S.compile(), "model_9.txt")

    def test_unit_event_test(self):
        A = BaseSpecies(1)
        A >> mobspy.Zero[1 / u.s]
        A(1 * u.mol)
        S = Simulation(A)
        S.level = -1
        with S.event_condition(A < 0.5 * u.mol):  # noqa: SIM300
            A(1 * u.mol)
        S.duration = 3
        assert compare_model(S.compile(), "model_10.txt")

    def test_event_all(self):
        Acka = BaseSpecies()
        Acka.a1, Acka.a2
        Baka = New(Acka)
        Baka.b1, Baka.b2
        Baka >> mobspy.Zero[1]
        S = Simulation(Baka)
        S.level = -1
        S.plot_data = False
        with S.event_time(0):
            set_counts({All[Baka]: 10, Baka.b1.a2: 20})
            Baka.a1.b1(30)
        assert compare_model(S.compile(), "model_24.txt")
        S.duration = 5
        S.step_size = 1
        S.run(plot_data=False)
        assert S.fres[Baka.a1.b1][0] == 30
        assert S.fres[Baka.b1.a2][0] == 20
        assert S.fres[Baka.a1.b2][0] == 10
        assert S.fres["Baka.a1.b1"][0] == 30
        assert S.fres["Baka.b1.a2"][0] == 20
        assert S.fres["Baka.a1.b2"][0] == 10
        for key in S.fres:
            if key == "Time":
                continue
            assert S.fres[key][-1] < 1

    def test_string_events_assignment(self):
        A = BaseSpecies()
        A.a1, A.a2, A.a3
        S = Simulation(A)
        S.level = -1
        with S.event_time(5):
            All[A](f"{A} + 1")
        with S.event_time(10):
            All[A.a1](f"{A} + 1")
        with S.event_time(15):
            A.a1(f"{A} + 1")
        compare_model(S.compile(), "model_30.txt")

    def test_event_reaction_not_allowed(self):
        try:
            A = BaseSpecies()
            A >> mobspy.Zero[1]
            S = Simulation(A)
            with S.event_time(0):
                mobspy.Zero >> A[1]
            assert False
        except (SystemExit, MobsPyError):
            assert True


class TestLogicOperators:
    def test_logic_operator_syntax(self):
        test_failed = False
        simlog.global_simlog_level = -1
        A, B = BaseSpecies(2)
        A.a1, A.a2, A.a3

        try:
            (10 >= A) >= 10  # noqa: SIM300
            test_failed = True
        except (SystemExit, MobsPyError):
            pass

        try:
            (10 >= A >= 10 >= A)
            test_failed = True
        except (SystemExit, MobsPyError):
            pass

        try:
            (10 >= A * A)  # noqa: SIM300
            test_failed = True
        except (SystemExit, MobsPyError):
            pass

        try:
            S1 = Simulation(A)
            S1.level = -1
            with S1.event_condition(B <= 10):
                A(100)
            S1.compile()
            test_failed = True
        except (SystemExit, MobsPyError):
            pass

        r1 = ((10 >= 2 * A) & (A <= 10)) | (10 >= A)  # noqa: SIM300
        S = Simulation(A)
        S.level = -1
        with S.event_condition(r1):
            A(100)
        assert compare_model(S.compile(), "model_16.txt")

        if test_failed:
            assert False

    def test_conditional_between_meta_species(self):
        Cu = BaseSpecies()
        Cu.c1, Cu.c2
        Azi, Byy = New(Cu)
        Azi.a1, Azi.a2, Byy.b1, Byy.b2
        Azi >> mobspy.Zero[1]
        Byy >> mobspy.Zero[0.1]
        Azi(200), Byy(50)
        S = Simulation(Azi | Byy)
        S.plot_data = False
        S.level = -1
        with S.event_condition(Azi.a1 <= Byy.b1):
            Azi(200)
        S.duration = 10
        S.method = "stochastic"
        assert compare_model(S.compile(), "model_19.txt")
        S.run(plot_data=False)
        for i in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
            assert S.fres[Azi][i] > S.fres[Byy][i]

    def test_conditional_between_meta_species_2(self):
        A, B = BaseSpecies()
        A >> mobspy.Zero[1]
        B >> mobspy.Zero[0.1]
        r1 = (A < B) & (A < B) | (A < B)
        A(200), B(50)
        S = Simulation(A | B)
        with S.event_condition(r1):
            A(200)
        S.level = -1
        assert compare_model(S.compile(), "model_20.txt")

    def test_bool_error(self):
        B = BaseSpecies()
        B >> mobspy.Zero[1]
        B(100)
        S = Simulation(B)
        simlog.global_simlog_level = -1
        try:
            with S.event_condition(B == 0):
                B(100)
            assert False
        except Exception:
            pass
        try:
            S.duration = True
            assert False
        except (SystemExit, MobsPyError):
            pass
        assert True
