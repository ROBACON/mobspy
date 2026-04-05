"""Tests for ODE syntax."""

from __future__ import annotations

from mobspy import BaseSpecies, New, Simulation, Zero
from mobspy.modules.functions import ms_exp
from mobspy.modules.ode_operator import dt


def test_ode_syntax_basic():
    A = BaseSpecies()
    dt[A] += -0.1 * A
    A(100)
    S = Simulation(A)
    S.compile(verbose=False)
    cm = S._concrete_model

    assert set(cm.species.keys()) == {"A"}
    assert cm.species["A"] == 100
    assert len(cm.reactions) == 1
    rxn = next(iter(cm.reactions.values()))
    assert rxn.reactants == [(1, "A")]
    assert rxn.products == [(2, "A")]
    assert rxn.kinetics == "(-0.1*A)"


def test_ode_syntax_two_species():
    A, B = BaseSpecies()
    dt[A] += -0.1 * A + 0.05 * B
    A(100), B(50)
    S = Simulation(A | B)
    S.compile(verbose=False)
    cm = S._concrete_model

    assert set(cm.species.keys()) == {"A", "B"}
    assert cm.species["A"] == 100
    assert cm.species["B"] == 50
    assert len(cm.reactions) == 1
    rxn = next(iter(cm.reactions.values()))
    assert sorted(rxn.reactants) == sorted([(1, "A"), (1, "B")])
    assert sorted(rxn.products) == sorted([(2, "A"), (1, "B")])
    assert rxn.kinetics == "((-0.1*A)+(0.05*B))"


def test_ode_applied_to_species():
    A = BaseSpecies()
    B = BaseSpecies()
    B.b1
    dt[A] += A
    dt[B.b1] += B.b1
    S = Simulation(A | B)
    S.compile(verbose=False)
    cm = S._concrete_model

    assert set(cm.species.keys()) == {"A", "B_dot_b1"}
    assert cm.species["A"] == 0
    assert cm.species["B_dot_b1"] == 0
    assert len(cm.reactions) == 2

    kinetics = {r.kinetics for r in cm.reactions.values()}
    assert "(1*A)" in kinetics
    assert "(1*B_dot_b1)" in kinetics

    for rxn in cm.reactions.values():
        if rxn.kinetics == "(1*A)":
            assert rxn.reactants == [(1, "A")]
            assert rxn.products == [(2, "A")]
        elif rxn.kinetics == "(1*B_dot_b1)":
            assert rxn.reactants == [(1, "B_dot_b1")]
            assert rxn.products == [(2, "B_dot_b1")]


def test_ode_neg_test():
    Neg, NegR = BaseSpecies()
    NegR.comp1
    dt[Neg] += -Neg
    dt[NegR] += -NegR.comp1
    S = Simulation(Neg | NegR)
    S.compile(verbose=False)
    cm = S._concrete_model

    assert set(cm.species.keys()) == {"Neg", "NegR_dot_comp1"}
    assert cm.species["Neg"] == 0
    assert cm.species["NegR_dot_comp1"] == 0
    assert len(cm.reactions) == 2

    kinetics = {r.kinetics for r in cm.reactions.values()}
    assert "(-1*Neg)" in kinetics
    assert "(-1*NegR_dot_comp1)" in kinetics


def test_ode_compartments():
    A = BaseSpecies()
    A.c1, A.c2
    Zero >> A.c1 @ 1
    A.c1 >> A.c2 @ 1
    dt[A] += -0.1 * A
    S = Simulation(A)
    S.compile(verbose=False)
    cm = S._concrete_model

    assert set(cm.species.keys()) == {"A_dot_c1", "A_dot_c2"}
    assert cm.species["A_dot_c1"] == 0
    assert cm.species["A_dot_c2"] == 0
    assert len(cm.reactions) == 4

    kinetics = sorted(r.kinetics for r in cm.reactions.values())
    assert "(-0.1*A_dot_c1)" in kinetics
    assert "(-0.1*A_dot_c2)" in kinetics

    # Creation reaction: Zero >> A.c1
    creation_rxns = [
        r
        for r in cm.reactions.values()
        if r.reactants == [] and r.products == [(1, "A_dot_c1")]
    ]
    assert len(creation_rxns) == 1
    assert "volume" in creation_rxns[0].kinetics

    # Transition reaction: A.c1 >> A.c2
    transition_rxns = [
        r
        for r in cm.reactions.values()
        if r.reactants == [(1, "A_dot_c1")] and r.products == [(1, "A_dot_c2")]
    ]
    assert len(transition_rxns) == 1


def test_ode_complex_expressions():
    A, B, C, D = BaseSpecies()
    dt[A] += 100 / (1 + B**2) - 0.1 * A
    dt[B] += (A * C) / (10 + A + C) - B / (5 + B)
    dt[C] += (A / (1 + A)) * (B / (1 + B)) - 0.05 * C * D
    dt[D] += (A**2 + B**2) / (100 + A**2 + B**2) * (1 - D / 1000)
    A(10), B(10), C(10), D(10)
    S = Simulation(A | B | C | D)
    S.compile(verbose=False)
    cm = S._concrete_model

    assert set(cm.species.keys()) == {"A", "B", "C", "D"}
    for sp in ("A", "B", "C", "D"):
        assert cm.species[sp] == 10
    assert len(cm.reactions) == 4

    all_kinetics = [r.kinetics for r in cm.reactions.values()]

    # Each ODE term produces one reaction; verify key subexpressions
    assert any("100" in k and "B^2" in k for k in all_kinetics)
    assert any("A*C" in k for k in all_kinetics)
    assert any("0.05*C" in k and "D" in k for k in all_kinetics)
    assert any("1000" in k and "D" in k for k in all_kinetics)


def test_ode_inheritance():
    Mortal = BaseSpecies()
    Human, Animal = New(Mortal)
    dt[Mortal] += -0.1 * Mortal
    Human(100), Animal(50)
    S = Simulation(Human | Animal)
    S.compile(verbose=False)
    cm = S._concrete_model

    assert set(cm.species.keys()) == {"Human", "Animal"}
    assert cm.species["Human"] == 100
    assert cm.species["Animal"] == 50
    assert len(cm.reactions) == 2

    kinetics = {r.kinetics for r in cm.reactions.values()}
    assert "(-0.1*Human)" in kinetics
    assert "(-0.1*Animal)" in kinetics


def test_ode_with_functions():
    A = BaseSpecies()
    dt[A] += 1 / (1 + ms_exp(A / 1000))
    A(100)
    S = Simulation(A)
    S.compile(verbose=False)
    cm = S._concrete_model

    assert set(cm.species.keys()) == {"A"}
    assert cm.species["A"] == 100
    assert len(cm.reactions) == 1
    rxn = next(iter(cm.reactions.values()))
    assert rxn.reactants == [(1, "A")]
    assert rxn.products == [(2, "A")]
    assert "exp" in rxn.kinetics
    assert "1000" in rxn.kinetics


def test_ode_neg():
    A = BaseSpecies()
    dt[A] -= A
    S = Simulation(A)
    S.compile(verbose=False)
    cm = S._concrete_model

    assert set(cm.species.keys()) == {"A"}
    assert cm.species["A"] == 0
    assert len(cm.reactions) == 1
    rxn = next(iter(cm.reactions.values()))
    # dt[A] -= A creates a reaction with stoichiometry 2 on reactant side
    assert rxn.reactants == [(2, "A")]
    assert rxn.products == [(1, "A")]
    assert rxn.kinetics == "(1*A)"
