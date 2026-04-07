"""Tests for species assignments (Assign context manager, .assign() method)."""

from __future__ import annotations

import mobspy
from mobspy import All, Assign, BaseSpecies, New, Simulation, logger, u
from mobspy.exceptions import MobsPyError


class TestAssign:
    def test_basic_assignment(self):
        A, B = BaseSpecies()
        A.assign(2 * B)
        B(100)
        S = Simulation(A | B)
        S.duration = 10
        S.step_size = 5
        S.level = -1
        S.run(plot_data=False)
        assert S.fres[A][-1] == 200

    def test_all_asgn_ops(self):
        A, B, C, D = BaseSpecies()
        A.a1.assign(B * 5)
        A.a2.assign(B + C)
        A.a3.assign(B - C)
        A.a4.assign(B / C)
        A.a5.assign(B / 200)
        A.a6.assign(B**2)
        A.a7.assign(B**D)

        B(200), C(100), D(2)
        S = Simulation(A | B | C | D)
        S.duration = 10
        S.step_size = 5
        S.level = -1
        S.run(plot_data=False)

        assert S.fres[A.a1][-1] == 1000
        assert S.fres[A.a2][-1] == 300
        assert S.fres[A.a3][-1] == 100
        assert S.fres[A.a4][-1] == 2
        assert S.fres[A.a5][-1] == 1
        assert S.fres[A.a6][-1] == 40000
        assert S.fres[A.a7][-1] == 40000

    def test_complex_assignments(self):
        A, B, C, D = BaseSpecies()
        A.assign(B * ((C + 5) / D))
        S = Simulation(A | B | C | D)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model
        assert cm is not None
        assert len(cm.species) == 4
        assert all(s in cm.species for s in ("A", "B", "C", "D"))
        assert len(cm.reactions) == 0
        assert len(cm.assignments) == 1
        asgn = next(iter(cm.assignments.values()))
        assert asgn.species == "A"
        for token in ("B", "C", "D", "5"):
            assert token in asgn.expression

    def test_assign_context_exit(self):
        try:
            logger.global_logger_level = -1
            A, B = BaseSpecies()
            A.assign(5 * B * (u.l / u.s) + 10 * B * (1 / u.s))
        except (SystemExit, MobsPyError):
            pass
        try:
            logger.global_logger_level = -1
            A >> mobspy.Zero @ 1
        except (SystemExit, MobsPyError):
            pass
        A, B = BaseSpecies()
        A >> mobspy.Zero @ 1
        B.assign(A / 2)
        S = Simulation(A | B)
        S.level = -1
        S.duration = 10
        S.step_size = 5
        S.run(plot_data=False)
        assert True

    def test_even_more_complex_assignments(self):
        Hi = BaseSpecies()
        Hi.h1, Hi.h2
        A, B, C, D = New(Hi)
        A.assign(B * ((C + 5) / D))
        S = Simulation(A | B | C | D)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model
        assert cm is not None
        # 4 base species * 2 characteristics (h1, h2) = 8 concrete species
        assert len(cm.species) == 8
        assert len(cm.reactions) == 0
        assert len(cm.assignments) == 1

    def test_assign_context_complex(self):
        A, B, C, D = BaseSpecies()
        B.b1, B.b2, B.b3
        C.c1, C.c2
        D.d1, D.d2
        with Assign:
            All[B]((C + D**2) * D)
        S = Simulation(A | B | C | D)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model
        assert cm is not None
        # A(1) + B.b1,B.b2,B.b3(3) + C.c1,C.c2(2) + D.d1,D.d2(2) = 8
        assert len(cm.species) == 8
        assert len(cm.reactions) == 0
        assert len(cm.assignments) == 3

    def test_assign_context_constant(self):
        A = BaseSpecies()
        with Assign:
            A(5)
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)
        cm = S._concrete_model
        assert cm is not None
        assert len(cm.species) == 1
        assert "A" in cm.species
        assert len(cm.reactions) == 0
        assert len(cm.assignments) == 1
        asgn = next(iter(cm.assignments.values()))
        assert "5" in asgn.expression
