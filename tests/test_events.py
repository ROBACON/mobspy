"""Tests for events (time-based, condition-based) and logic operators."""

from __future__ import annotations

import mobspy
from mobspy import All, BaseSpecies, New, Simulation, logger, set_counts, u
from mobspy.exceptions import MobsPyError


class TestTimeEvents:
    def test_event_type(self):
        A, B, C, D, E, F = BaseSpecies(6)
        A + B >> mobspy.Zero @ 1
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
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 6
        assert cm.species["A"] == 50
        assert cm.species["B"] == 50
        assert cm.species["C"] == 0
        assert cm.species["F"] == 0

        real_rxns = {
            k: v for k, v in cm.reactions.items() if not k.startswith("phantom_")
        }
        assert len(real_rxns) == 1
        rxn = next(iter(real_rxns.values()))
        assert ("A" in rxn.kinetics) and ("B" in rxn.kinetics)

        assert len(cm.events) == 4
        triggers = [e.trigger for e in cm.events.values()]
        assignments_targets = [t for e in cm.events.values() for t, _ in e.assignments]
        assert any("<=" in t and "A" in t and "B" in t for t in triggers)
        assert any("B" in t and "<=" in t and "A" not in t for t in triggers)
        assert any(t == "true" for t in triggers)
        assert "C" in assignments_targets
        assert "D" in assignments_targets
        assert "E" in assignments_targets
        assert "F" in assignments_targets

    def test_reacting_species_event(self):
        B = BaseSpecies(1)
        B.b1, B.b2
        A = New(B)
        B >> mobspy.Zero @ 1
        A.b2 >> mobspy.Zero @ 0.5
        A.a1 >> mobspy.Zero @ 1
        A.a1(100), A.b2(100), B.b1(100)
        S = Simulation(A | B)
        S.level = -1
        with S.event_condition((A.a1 <= 10) & (B.b1 <= 10)):
            A.a1(100)
        S.duration = 5
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 4
        assert "A.a1.b1" in cm.species or "A.b1.a1" in cm.species
        assert "B.b1" in cm.species
        assert "B.b2" in cm.species

        assert len(cm.reactions) == 7

        assert len(cm.events) == 1
        event = next(iter(cm.events.values()))
        assert "<=" in event.trigger
        assert "10" in event.trigger
        assert len(event.assignments) == 1
        target, value = event.assignments[0]
        assert "A" in target
        assert str(value) == "100"

    def test_unit_event_test(self):
        A = BaseSpecies(1)
        A >> mobspy.Zero @ (1 / u.s)
        A(1 * u.mol)
        S = Simulation(A)
        S.level = -1
        with S.event_condition(A < 0.5 * u.mol):  # noqa: SIM300
            A(1 * u.mol)
        S.duration = 3
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 1
        assert "A" in cm.species

        assert len(cm.reactions) == 1

        assert len(cm.events) == 1
        event = next(iter(cm.events.values()))
        assert "<" in event.trigger
        assert "0.5" in event.trigger or "A" in event.trigger
        assert len(event.assignments) == 1
        target, value = event.assignments[0]
        assert target == "A"

    def test_event_all(self):
        Acka = BaseSpecies()
        Acka.a1, Acka.a2
        Baka = New(Acka)
        Baka.b1, Baka.b2
        Baka >> mobspy.Zero @ 1
        S = Simulation(Baka)
        S.level = -1
        S.plot_data = False
        with S.event_time(0):
            set_counts({All[Baka]: 10, Baka.b1.a2: 20})
            Baka.a1.b1(30)
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 4
        assert "Baka.a1.b1" in cm.species
        assert "Baka.a1.b2" in cm.species
        assert "Baka.a2.b1" in cm.species
        assert "Baka.a2.b2" in cm.species

        assert len(cm.reactions) == 4

        assert len(cm.events) == 1
        event = next(iter(cm.events.values()))
        assert event.trigger == "true"
        assert len(event.assignments) == 4
        assignment_dict = dict(event.assignments)
        assert str(assignment_dict["Baka_dot_a1_dot_b1"]) == "30"
        assert str(assignment_dict["Baka_dot_a2_dot_b1"]) == "20"
        assert str(assignment_dict["Baka_dot_a1_dot_b2"]) == "10"
        assert str(assignment_dict["Baka_dot_a2_dot_b2"]) == "10"

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
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 3
        assert "A.a1" in cm.species
        assert "A.a2" in cm.species
        assert "A.a3" in cm.species

        real_rxns = {
            k: v for k, v in cm.reactions.items() if not k.startswith("phantom_")
        }
        assert len(real_rxns) == 0
        assert len(cm.events) == 3

        delays = sorted(str(e.delay) for e in cm.events.values())
        assert "10" in delays
        assert "15" in delays
        assert "5" in delays

        for event in cm.events.values():
            assert event.trigger == "true"
            for _, expr in event.assignments:
                assert "+ 1" in str(expr)

        # The event at delay=5 assigns all three species
        event_5 = next(e for e in cm.events.values() if str(e.delay) == "5")
        assert len(event_5.assignments) == 3

        # The event at delay=10 assigns only A.a1 variants
        event_10 = next(e for e in cm.events.values() if str(e.delay) == "10")
        assert len(event_10.assignments) == 1
        assert event_10.assignments[0][0] == "A_dot_a1"

        # The event at delay=15 assigns only A.a1
        event_15 = next(e for e in cm.events.values() if str(e.delay) == "15")
        assert len(event_15.assignments) == 1
        assert event_15.assignments[0][0] == "A_dot_a1"

    def test_event_reaction_not_allowed(self):
        try:
            A = BaseSpecies()
            A >> mobspy.Zero @ 1
            S = Simulation(A)
            with S.event_time(0):
                mobspy.Zero >> A @ 1
            assert False
        except (SystemExit, MobsPyError):
            assert True


class TestLogicOperators:
    def test_logic_operator_syntax(self):
        test_failed = False
        logger.global_logger_level = -1
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
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 3
        assert "A.a1" in cm.species
        assert "A.a2" in cm.species
        assert "A.a3" in cm.species

        real_rxns = {
            k: v for k, v in cm.reactions.items() if not k.startswith("phantom_")
        }
        assert len(real_rxns) == 0
        assert len(cm.events) == 1

        event = next(iter(cm.events.values()))
        assert "<=" in event.trigger or ">=" in event.trigger
        assert "&&" in event.trigger or "||" in event.trigger
        assert len(event.assignments) == 1
        assert str(event.assignments[0][1]) == "100"

        if test_failed:
            assert False

    def test_conditional_between_meta_species(self):
        Cu = BaseSpecies()
        Cu.c1, Cu.c2
        Azi, Byy = New(Cu)
        Azi.a1, Azi.a2, Byy.b1, Byy.b2
        Azi >> mobspy.Zero @ 1
        Byy >> mobspy.Zero @ 0.1
        Azi(200), Byy(50)
        S = Simulation(Azi | Byy)
        S.plot_data = False
        S.level = -1
        with S.event_condition(Azi.a1 <= Byy.b1):
            Azi(200)
        S.duration = 10
        S.method = "stochastic"
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 8
        assert sum(1 for k in cm.species if k.startswith("Azi")) == 4
        assert sum(1 for k in cm.species if k.startswith("Byy")) == 4

        assert len(cm.reactions) == 8
        azi_rxns = [r for r in cm.reactions.values() if "Azi" in r.kinetics]
        byy_rxns = [r for r in cm.reactions.values() if "Byy" in r.kinetics]
        assert len(azi_rxns) == 4
        assert len(byy_rxns) == 4
        assert all("* 1" in r.kinetics for r in azi_rxns)
        assert all("* 0.1" in r.kinetics for r in byy_rxns)

        assert len(cm.events) == 1
        event = next(iter(cm.events.values()))
        assert "<=" in event.trigger
        assert "Azi" in event.trigger
        assert "Byy" in event.trigger

        S.run(plot_data=False)
        for i in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
            assert S.fres[Azi][i] > S.fres[Byy][i]

    def test_conditional_between_meta_species_2(self):
        A, B = BaseSpecies()
        A >> mobspy.Zero @ 1
        B >> mobspy.Zero @ 0.1
        r1 = (A < B) & (A < B) | (A < B)
        A(200), B(50)
        S = Simulation(A | B)
        with S.event_condition(r1):
            A(200)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 2
        assert cm.species["A"] == 200
        assert cm.species["B"] == 50

        assert len(cm.reactions) == 2

        assert len(cm.events) == 1
        event = next(iter(cm.events.values()))
        assert "<" in event.trigger
        assert "A" in event.trigger
        assert "B" in event.trigger
        assert "&&" in event.trigger or "||" in event.trigger
        assert len(event.assignments) == 1
        assert event.assignments[0][0] == "A"
        assert str(event.assignments[0][1]) == "200"

    def test_bool_error(self):
        B = BaseSpecies()
        B >> mobspy.Zero @ 1
        B(100)
        S = Simulation(B)
        logger.global_logger_level = -1
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
