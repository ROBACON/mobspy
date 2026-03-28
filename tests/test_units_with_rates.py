"""Tests for units combined with non-mass-action rates, is_a queries,
characteristic dot notation, math functions, and Count/Concentration operators.

These tests exercise the full pipeline: expression AST construction,
unit resolution, and SBML string generation.
"""

from __future__ import annotations

import re

import pytest

import mobspy
from mobspy import BaseSpecies, New, Simulation, u
from mobspy.modules.functions import ms_abs, ms_exp, ms_logn, ms_sqrt
from mobspy.modules.ode_operator import dt

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _kinetics_from_compiled(compiled: str) -> dict[str, str]:
    """Extract reaction_name -> kinetics string from compiled output."""
    result = {}
    for line in compiled.splitlines():
        m = re.match(r"(reaction_\d+),\{.*'kin':\s*'([^']+)'", line)
        if m:
            result[m.group(1)] = m.group(2)
    return result


def _species_counts_from_compiled(compiled: str) -> dict[str, str]:
    """Extract species -> initial count from compiled output."""
    result = {}
    in_species = False
    for line in compiled.splitlines():
        if line.strip() == "Species":
            in_species = True
            continue
        if in_species and line.strip() == "":
            break
        if in_species and "," in line:
            name, count = line.strip().split(",", 1)
            result[name] = count
    return result


# ===========================================================================
# Unit-aware rate functions
# ===========================================================================


class TestUnitRateFunctions:
    """Rate functions that return Pint Quantity values."""

    def test_first_order_rate_with_time_unit(self):
        A = BaseSpecies()
        A >> mobspy.Zero[lambda r: 2 / u.hour * r]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 1
        rate_str = next(iter(kin.values()))
        assert "A" in rate_str
        assert "0.000555" in rate_str or "5.555" in rate_str

    def test_second_order_concentration_rate(self):
        A, B = BaseSpecies()
        (
            A + B
            >> mobspy.Zero[
                lambda r1, r2: (
                    (1 * u.millimolar / u.hour)
                    * (1 + 10 * u.millimolar / r1 + 20 * u.millimolar / r2)
                )
            ]
        )
        A(100), B(200)
        S = Simulation(A | B)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 1
        rate_str = next(iter(kin.values()))
        assert "volume" in rate_str

    def test_zero_order_unit_rate(self):
        A = BaseSpecies()
        mobspy.Zero >> A[1 / u.s]
        A(0)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        rate_str = next(iter(kin.values()))
        assert "1.0" in rate_str or "1" in rate_str

    def test_unit_rate_with_volume(self):
        A = BaseSpecies()
        A + A >> mobspy.Zero[1 * u.liter / u.second]
        A(100)
        S = Simulation(A)
        S.level = -1
        S.volume = 1 * u.liter
        result = S.compile()
        assert "volume" in result

    def test_rate_function_division_by_species(self):
        A, B = BaseSpecies()
        A >> mobspy.Zero[lambda r: (10 / u.s) / (1 + r)]
        A(50)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        rate_str = next(iter(kin.values()))
        assert "A" in rate_str

    def test_mixed_unit_and_unitless_in_expression(self):
        A = BaseSpecies()
        A >> mobspy.Zero[lambda r: (1 / u.minute) * (r / (50 + r))]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 1

    def test_rate_function_minute_hour_conversion(self):
        """Verify that 60/u.minute and 1/u.second produce the same rate."""
        A = BaseSpecies()
        A >> mobspy.Zero[60 / u.minute]
        A(100)
        S1 = Simulation(A)
        S1.level = -1
        r1 = S1.compile()

        A.reset_reactions()
        A >> mobspy.Zero[1 / u.second]
        S2 = Simulation(A)
        S2.level = -1
        r2 = S2.compile()

        kin1 = _kinetics_from_compiled(r1)
        kin2 = _kinetics_from_compiled(r2)
        assert kin1[next(iter(kin1.keys()))] == kin2[next(iter(kin2.keys()))]


# ===========================================================================
# is_a queries in rate functions
# ===========================================================================


class TestIsAInRates:
    """Test is_a() queries inside rate function lambdas."""

    def test_is_a_basic(self):
        Mortal = BaseSpecies()
        Protein = New(Mortal)
        mRNA = New(Mortal)
        Mortal >> mobspy.Zero[lambda r: 0.01 if r.is_a(Protein) else 1.0]
        Protein(100), mRNA(50)
        S = Simulation(Protein | mRNA)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        rates = list(kin.values())
        assert len(rates) == 2
        rate_values = sorted(rates)
        assert "0.01" in rate_values[0]
        assert "1.0" in rate_values[1] or "1" in rate_values[1]

    def test_is_a_with_units(self):
        Mortal = BaseSpecies()
        Fast = New(Mortal)
        Slow = New(Mortal)
        Mortal >> mobspy.Zero[lambda r: 10 / u.hour if r.is_a(Fast) else 1 / u.hour]
        Fast(100), Slow(100)
        S = Simulation(Fast | Slow)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        rates = list(kin.values())
        assert len(rates) == 2
        # Both should have been converted from /hour to /second
        for r in rates:
            assert "volume" not in r  # first order, no volume needed

    def test_is_a_in_second_order_with_units(self):
        Base = BaseSpecies()
        TypeA = New(Base)
        TypeB = New(Base)
        (
            Base + Base
            >> mobspy.Zero[
                lambda r1, r2: 1 * u.l / u.s if r1.is_a(TypeA) else 0.001 * u.l / u.s
            ]
        )
        TypeA(100), TypeB(50)
        S = Simulation(TypeA | TypeB)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        # TypeA+TypeA and TypeA+TypeB get rate 1, TypeB+TypeB and TypeB+TypeA get 0.001
        assert len(kin) >= 2

    def test_is_a_combined_with_arithmetic(self):
        Organism = BaseSpecies()
        Bacterium = New(Organism)
        Virus = New(Organism)
        (
            Organism
            >> mobspy.Zero[
                lambda r: (0.5 / u.min) * r if r.is_a(Bacterium) else (2 / u.min) * r
            ]
        )
        Bacterium(1000), Virus(500)
        S = Simulation(Bacterium | Virus)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 2

    def test_is_a_three_levels(self):
        Base = BaseSpecies()
        Mid = New(Base)
        Leaf = New(Mid)
        # Leaf.is_a(Leaf) -> True, Leaf.is_a(Mid) -> True
        # Mid.is_a(Leaf) -> False, Mid.is_a(Mid) -> True
        # Base has no concrete instances (not in sim) so only Mid and Leaf get reactions
        Base >> mobspy.Zero[lambda r: 1 if r.is_a(Leaf) else 0.5]
        Leaf(10), Mid(20)
        S = Simulation(Leaf | Mid)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        # Leaf gets rate 1 (mass-action: Leaf * 1), Mid gets rate 0.5 (Mid * 0.5)
        rates = sorted(kin.values())
        assert len(rates) == 2
        assert any("Leaf" in r and "1" in r for r in rates)
        assert any("Mid" in r and "0.5" in r for r in rates)


# ===========================================================================
# Characteristic dot notation in rate functions
# ===========================================================================


class TestCharacteristicQueries:
    """Test dot notation (r.alive, r.state) in rate lambdas."""

    def test_single_characteristic(self):
        A = BaseSpecies()
        A.alive, A.dead
        A >> mobspy.Zero[lambda r: 1 if r.alive else 0.5]
        A.alive(50), A.dead(50)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        # Mass-action rates: alive -> "A.alive * 1", dead -> "A.dead * 0.5"
        rates = sorted(kin.values())
        assert len(rates) == 2
        assert any("A.dead" in r and "0.5" in r for r in rates)
        assert any("A.alive" in r and "1" in r for r in rates)

    def test_characteristic_with_units(self):
        Cell = BaseSpecies()
        Cell.healthy, Cell.infected
        Cell >> mobspy.Zero[lambda r: 0.01 / u.hour if r.healthy else 1 / u.hour]
        Cell.healthy(1000), Cell.infected(100)
        S = Simulation(Cell)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 2

    def test_two_characteristics_product(self):
        A, B = BaseSpecies()
        A.a1, A.a2
        B.b1, B.b2
        Combo = A * B
        Combo >> mobspy.Zero[lambda r: 1 if r.a1 else (2 if r.b1 else 3)]
        S = Simulation(Combo)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 4  # a1.b1, a1.b2, a2.b1, a2.b2

    def test_characteristic_with_rate_arithmetic_and_units(self):
        Gene = BaseSpecies()
        Gene.on, Gene.off
        Gene >> mobspy.Zero[lambda r: (5 / u.min) * r if r.on else (0.1 / u.min) * r]
        Gene.on(10), Gene.off(90)
        S = Simulation(Gene)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 2
        for rate_str in kin.values():
            assert "Gene" in rate_str

    def test_characteristic_in_two_reactant_rate(self):
        A, B = BaseSpecies()
        A.fast, A.slow
        A + B >> mobspy.Zero[lambda r1, r2: 1 if r1.fast else 0.1]
        A.fast(50), A.slow(50), B(100)
        S = Simulation(A | B)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 2

    def test_is_a_and_characteristic_combined(self):
        Mortal = BaseSpecies()
        Mortal.young, Mortal.old
        Human = New(Mortal)
        Animal = New(Mortal)
        Mortal >> mobspy.Zero[lambda r: 0.1 if r.is_a(Human) and r.young else 1.0]
        Human.young(50), Human.old(50)
        Animal.young(30), Animal.old(30)
        S = Simulation(Human | Animal)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        # Human.young -> 0.1 (mass-action: "Human.young * 0.1")
        # Human.old, Animal.young, Animal.old -> 1.0 (mass-action: "X * 1")
        assert len(kin) == 4
        rates = sorted(kin.values())
        assert any("0.1" in r for r in rates)
        assert sum(1 for r in rates if "* 1" in r) == 3

    def test_is_a_and_characteristic_with_units(self):
        Base = BaseSpecies()
        Base.active, Base.inactive
        TypeA = New(Base)
        TypeB = New(Base)
        (
            Base
            >> mobspy.Zero[
                lambda r: 10 / u.s if r.is_a(TypeA) and r.active else 0.1 / u.s
            ]
        )
        TypeA.active(100), TypeA.inactive(50)
        TypeB.active(80), TypeB.inactive(40)
        S = Simulation(TypeA | TypeB)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        # TypeA.active, TypeA.inactive, TypeB.active, TypeB.inactive
        assert len(kin) == 4


# ===========================================================================
# ODE syntax with units
# ===========================================================================


class TestODEWithUnits:
    """ODE dt[] syntax combined with units and math functions."""

    def test_ode_basic_species_arithmetic(self):
        A, B = BaseSpecies()
        dt[A] += -0.1 * A + 0.05 * B
        dt[B] += 0.1 * A - 0.05 * B
        A(100), B(0)
        S = Simulation(A | B)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) >= 2

    def test_ode_with_exp_function(self):
        A = BaseSpecies()
        dt[A] += 10 / (1 + ms_exp(A / 100)) - 0.1 * A
        A(50)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        rate_strs = list(kin.values())
        assert any("exp" in r for r in rate_strs)

    def test_ode_with_log_function(self):
        A = BaseSpecies()
        dt[A] += ms_logn(1 + A) - 0.01 * A
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert any("log" in r for r in kin.values())

    def test_ode_with_sqrt_function(self):
        A = BaseSpecies()
        dt[A] += ms_sqrt(A) - 0.1 * A
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert any("sqrt" in r for r in kin.values())

    def test_ode_with_abs_function(self):
        A, B = BaseSpecies()
        dt[A] += ms_abs(A - B) - 0.05 * A
        A(100), B(50)
        S = Simulation(A | B)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert any("abs" in r for r in kin.values())

    def test_ode_death_syntax(self):
        A = BaseSpecies()
        dt[A] -= 0.5 * A
        A(200)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 1

    def test_ode_hill_function(self):
        A, B = BaseSpecies()
        n = 4
        K = 100
        dt[A] += B**n / (K**n + B**n) - 0.1 * A
        A(10), B(50)
        S = Simulation(A | B)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert any("^" in r for r in kin.values())

    def test_ode_multiple_species_complex(self):
        A, B, C = BaseSpecies()
        dt[A] += (B * C) / (1 + B + C) - 0.1 * A
        dt[B] += A / (1 + A) - B / (10 + B)
        dt[C] += ms_exp(-A / 100) * B - 0.05 * C
        A(10), B(20), C(30)
        S = Simulation(A | B | C)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) >= 3

    def test_ode_with_inheritance(self):
        Mortal = BaseSpecies()
        Human = New(Mortal)
        Animal = New(Mortal)
        dt[Mortal] += -0.1 * Mortal
        Human(100), Animal(50)
        S = Simulation(Human | Animal)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) >= 2

    def test_ode_with_characteristics(self):
        A = BaseSpecies()
        A.c1, A.c2
        dt[A] += -0.1 * A
        A.c1(100), A.c2(50)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) >= 2


# ===========================================================================
# Parameters with units
# ===========================================================================


class TestParametersWithUnits:
    """ModelParameters combined with units in rate expressions."""

    def test_parameter_with_hour_unit(self):
        from mobspy import ModelParameters

        A = BaseSpecies()
        k = ModelParameters(1 / u.h)
        A >> mobspy.Zero[k]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert "k" in result

    def test_parameter_in_expression_with_units(self):
        from mobspy import ModelParameters

        A = BaseSpecies()
        k1, k2 = ModelParameters(0.1 / u.s, 50)
        A >> mobspy.Zero[lambda r: k1 * r / (k2 + r)]
        A(200)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        rate_str = next(iter(kin.values()))
        assert "k1" in rate_str
        assert "k2" in rate_str
        assert "A" in rate_str

    def test_parameter_concentration_unit(self):
        from mobspy import ModelParameters

        A = BaseSpecies()
        Km = ModelParameters(10 * u.millimolar)
        Vmax = ModelParameters(5 / u.s)
        A >> mobspy.Zero[lambda r: Vmax * r / (Km + r)]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 1


# ===========================================================================
# Unit counts and concentrations
# ===========================================================================


class TestUnitCounts:
    """Verify initial count conversion with various unit formats."""

    def test_count_in_molecules(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1]
        A(1000)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        counts = _species_counts_from_compiled(result)
        assert counts["A"] == "1000"

    def test_count_in_moles(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1]
        A(1 * u.mol)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        counts = _species_counts_from_compiled(result)
        count_val = float(counts["A"])
        # With model_context, species in moles stay as moles (1.0)
        assert count_val == pytest.approx(1.0, rel=1e-5)

    def test_count_in_concentration_with_volume(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1]
        A(1 * u.mol / u.liter)
        S = Simulation(A)
        S.level = -1
        S.volume = 1 * u.liter
        result = S.compile()
        counts = _species_counts_from_compiled(result)
        count_val = float(counts["A"])
        # With model_context, 1 mol/L * 1 L = 1 mol (stays in moles)
        assert count_val == pytest.approx(1.0, rel=1e-5)

    def test_dimensionless_count(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1]
        a = 100 * u.l / u.l
        A(a)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        counts = _species_counts_from_compiled(result)
        assert float(counts["A"]) == pytest.approx(100.0)


# ===========================================================================
# Volume and dimension consistency
# ===========================================================================


class TestVolumeAndDimension:
    """Tests for volume specification and dimensional consistency."""

    def test_2d_model_with_units(self):
        A = BaseSpecies()
        A + A >> mobspy.Zero[1 * u.decimeter**2 / u.s]
        A(100)
        S = Simulation(A)
        S.level = -1
        S.volume = 1 * u.m**2
        result = S.compile()
        assert "volume" in result

    def test_3d_model_default(self):
        A = BaseSpecies()
        A + A >> mobspy.Zero[1 * u.liter / u.second]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert "volume" in result

    def test_volume_conversion_ml_to_l(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1]
        A(100)
        S = Simulation(A)
        S.volume = 500 * u.mL
        S.level = -1
        result = S.compile()
        # 500 mL = 0.5 L = 0.5 dm^3
        assert "0.5" in result


# ===========================================================================
# Reversible reactions with units
# ===========================================================================


class TestReversibleWithUnits:
    """Reversible reactions with unit-aware rates."""

    def test_rev_constant_rates_with_units(self):
        from mobspy import Rev

        A, B = BaseSpecies()
        Rev[A >> B][1 / u.s, 2 / u.s]
        A(100), B(0)
        S = Simulation(A | B)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 2

    def test_rev_lambda_rates_with_units(self):
        from mobspy import Rev

        A, B = BaseSpecies()
        Rev[A >> B][
            lambda r: (1 / u.hour) * r,
            lambda r: (0.5 / u.hour) * r,
        ]
        A(100), B(0)
        S = Simulation(A | B)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 2


# ===========================================================================
# 2D rate queries with characteristics
# ===========================================================================


class TestCharacteristicEqualityInRates:
    """Tests using characteristic equality checks like Location(r1) == Location(r2)."""

    def test_location_equality_with_units(self):
        Color, Location = BaseSpecies()
        Color.red, Color.blue
        Location.here, Location.there
        Something = Color * Location
        rate = lambda r1, r2: (
            1 * u.decimeter**2 / u.h
            if Location(r1) == Location(r2)
            else 0.5 * u.decimeter**2 / u.h
        )
        2 * Something >> 3 * Something[rate]
        S = Simulation(Something)
        S.level = -1
        S.volume = 1 * u.m**2
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) > 0


# ===========================================================================
# Edge cases and error handling
# ===========================================================================


class TestEdgeCases:
    """Edge cases in unit handling."""

    def test_zero_rate_with_units(self):
        """Zero rate with units should compile without error."""
        A = BaseSpecies()
        A >> mobspy.Zero[0 / u.s]
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        # Zero-rate reactions may be omitted or have rate "0"
        kin = _kinetics_from_compiled(result)
        if len(kin) > 0:
            assert kin[next(iter(kin.keys()))] == "0"

    def test_very_small_rate_with_units(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1e-15 / u.s]
        A(1000)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 1

    def test_rate_function_returning_species_times_unit(self):
        A = BaseSpecies()
        A >> 2 * A[lambda r: r * (1 / u.s)]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        kin = _kinetics_from_compiled(result)
        assert len(kin) == 1

    def test_multiple_unit_systems_consistent(self):
        """Same rate expressed in different units should compile."""
        A = BaseSpecies()
        A >> mobspy.Zero[3600 / u.hour]
        A(100)
        S1 = Simulation(A)
        S1.level = -1
        r1 = S1.compile()

        A.reset_reactions()
        A >> mobspy.Zero[60 / u.minute]
        S2 = Simulation(A)
        S2.level = -1
        r2 = S2.compile()

        kin1 = next(iter(_kinetics_from_compiled(r1).values()))
        kin2 = next(iter(_kinetics_from_compiled(r2).values()))
        assert kin1 == kin2

    def test_event_with_unit_condition(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1 / u.s]
        A(1 * u.mol)
        S = Simulation(A)
        S.level = -1
        with S.event_condition(A < 0.5 * u.mol):  # noqa: SIM300
            A(1 * u.mol)
        S.duration = 3
        result = S.compile()
        assert "Events" in result or "event" in result.lower()


# ===========================================================================
# Simulation run validation (numerical correctness)
# ===========================================================================


class TestSimulationRunWithUnits:
    """Actually run simulations and check numerical output."""

    def test_exponential_decay_with_units(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1 / u.s]
        A(1000)
        S = Simulation(A)
        S.duration = 5
        S.level = -1
        S.run(plot_data=False)
        # After 5 seconds at rate 1/s, should be ~1000*exp(-5) ~ 6.7
        final = S.fres[A][-1]
        assert final < 50

    def test_growth_with_unit_rate(self):
        A = BaseSpecies()
        A >> 2 * A[1 / u.minute]
        A(10)
        S = Simulation(A)
        S.duration = 60  # 60 seconds = 1 minute
        S.level = -1
        S.run(plot_data=False)
        final = S.fres[A][-1]
        # Should have grown: 10 * exp(1/60 * 60) ~ 10*e ~ 27
        assert final > 20

    def test_is_a_produces_different_dynamics(self):
        Mortal = BaseSpecies()
        Fast = New(Mortal)
        Slow = New(Mortal)
        Mortal >> mobspy.Zero[lambda r: 10 / u.s if r.is_a(Fast) else 0.001 / u.s]
        Fast(1000), Slow(1000)
        S = Simulation(Fast | Slow)
        S.duration = 1
        S.level = -1
        S.run(plot_data=False)
        # Fast should decay much more than Slow
        assert S.fres[Fast][-1] < S.fres[Slow][-1]

    def test_characteristic_rate_produces_different_dynamics(self):
        Cell = BaseSpecies()
        Cell.healthy, Cell.sick
        Cell >> mobspy.Zero[lambda r: 0.001 if r.healthy else 1]
        Cell.healthy(500), Cell.sick(500)
        S = Simulation(Cell)
        S.duration = 5
        S.level = -1
        S.run(plot_data=False)
        # Sick cells decay faster
        assert S.fres["Cell.healthy"][-1] > S.fres["Cell.sick"][-1]


# ===========================================================================
# End-to-end: model runs in user-provided units
# ===========================================================================


class TestModelInUserUnits:
    """Verify that models with substance units run in the user's unit system."""

    def test_compile_molar_species_uses_moles(self):
        """Species in moles should stay in moles, not be converted to N_A counts."""
        A = BaseSpecies()
        A >> mobspy.Zero[1 / u.minute]
        A(2 * u.mmol)
        S = Simulation(A)
        S.level = -1
        compiled = S.compile()
        counts = _species_counts_from_compiled(compiled)
        # 2 mmol stays as 2 (model substance unit = millimole, first-unit-wins)
        assert abs(float(counts["A"]) - 2.0) < 1e-6

    def test_compile_molar_concentration_with_volume(self):
        """Molar concentration * volume should give moles."""
        A = BaseSpecies()
        A >> mobspy.Zero[1 / u.second]
        A(0.5 * u.mol / u.liter)
        S = Simulation(A)
        S.volume = 2 * u.liter
        S.level = -1
        compiled = S.compile()
        counts = _species_counts_from_compiled(compiled)
        # 0.5 M * 2 L = 1.0 mol
        assert abs(float(counts["A"]) - 1.0) < 1e-6

    def test_run_molar_model_time_axis(self):
        """Time axis in model time units when duration uses minutes."""
        A = BaseSpecies()
        A >> mobspy.Zero[0.1 / u.minute]
        A(1 * u.mol)
        S = Simulation(A)
        S.duration = 10 * u.minute
        S.level = -1
        S.run(plot_data=False)
        time_vals = S.fres["Time"]
        # Time should go from 0 to 10 (minutes), not 0 to 600 (seconds)
        assert time_vals[-1] == pytest.approx(10.0, rel=1e-2)

    def test_run_molar_model_dynamics(self):
        """Exponential decay A -> 0 with rate k: A(t) = A0 * exp(-k*t)."""
        A = BaseSpecies()
        k = 0.5  # 1/second (dimensionless)
        A >> mobspy.Zero[k]
        A(1 * u.mol)
        S = Simulation(A)
        S.duration = 4
        S.level = -1
        S.run(plot_data=False)
        final_A = S.fres[A][-1]
        import math

        expected = 1.0 * math.exp(-k * 4)
        assert final_A == pytest.approx(expected, rel=0.05)

    def test_plain_number_model_unchanged(self):
        """Model with no Pint units should produce identical results to legacy."""
        A = BaseSpecies()
        A >> mobspy.Zero[0.1]
        A(100)
        S = Simulation(A)
        S.duration = 10
        S.level = -1
        S.run(plot_data=False)
        final_A = S.fres[A][-1]
        import math

        expected = 100 * math.exp(-0.1 * 10)
        assert final_A == pytest.approx(expected, rel=0.05)
