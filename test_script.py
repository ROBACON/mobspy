# This script tests the model creation and compilation
# It also test the calculation capabilities
# It uses some simple model and assertions

# import pytest
from mobspy import *
import sys
import matplotlib.pyplot as plt


# TODO Plot has random order for species names


# Compare results with expected file
def compare_model(comp_results, file_name):
    with open(file_name, 'r') as file:
        for r, line in zip(comp_results.split('\n'), file.readlines()):
            line = line.replace('\n', '')
            if r != line:
                return False
    return True


# Model to test the basics
def test_model_1():
    A, B, C = BaseSpecies(3)
    A + B >> C[1]
    MySim = Simulation(A | B | C)
    MySim.level = 0
    results = MySim.compile()
    assert compare_model(results, 'test_tools/model_1.txt')


# Model to test basic inheritance
def test_model_2():
    Carnivore, Herbivore = BaseSpecies(2)
    Cat, Dog = New(Carnivore, 2)
    Carnivore + Herbivore(1 * u.mol) >> Carnivore[1]
    Cat(1 * u.mol), Dog(1 * u.mol)
    MySim = Simulation(Cat | Dog | Herbivore)
    MySim.level = 0
    MySim.volume = 1 * u.meter ** 2
    results = MySim.compile()
    assert compare_model(results, 'test_tools/model_2.txt')


# Model to test species multiplication
def test_model_3():
    MGMT, Blue_Oyster_Cult, The_Smiths = BaseSpecies(3)
    MGMT.eletric_fell, MGMT.little_dark_age, MGMT.kids
    Blue_Oyster_Cult.burning_for_you >> Blue_Oyster_Cult.reaper[1]
    The_Smiths.stop_me >> The_Smiths.charming_man[1]
    Music = MGMT * Blue_Oyster_Cult * The_Smiths
    MySim = Simulation(Music)
    MySim.level = 0
    results = MySim.compile()
    assert compare_model(results, 'test_tools/model_3.txt')


# Model to test inheritance queries
# All bacterias are infected by any virus here
def test_model_4():
    Bacteria, Virus = BaseSpecies(2)
    B1, B2 = New(Bacteria, 2)
    V1, V2 = New(Virus, 2)
    Bacteria.not_infected + Virus >> Bacteria.infected[1]
    MySim = Simulation(B1 | B2 | V1 | V2)
    MySim.level = 0
    results = MySim.compile()
    assert compare_model(results, 'test_tools/model_4.txt')


# Model to test round-robin and stoichiometry
def test_model_5():
    A = BaseSpecies(1)
    B, C = New(A, 2)
    A >> 2 * A[1]
    2 * A >> 3 * A[1]
    MySim = Simulation(B | C)
    MySim.level = 0
    results = MySim.compile()
    assert compare_model(results, 'test_tools/model_5.txt')


def test_model_6():
    # This model tests species that are not referenced in the reactants (we call them Born Species)
    A = BaseSpecies(1)
    B = New(A)
    C = New(A)
    B.b1, B.b2, C.c1, C.c2
    Zero >> 2 * A[1]
    MySim = Simulation(B | C)
    results = MySim.compile()
    assert compare_model(results, 'test_tools/model_6.txt')


def test_model_7():
    def oscillator(beta_m=5, beta_p=10, gamma_m=1, gamma_p=0.01, k=1, n=4, leaky=0.0001):
        Mortal, Creator = BaseSpecies(2)

        mRNA = Mortal * Creator
        Protein = New(Mortal)

        # Repression reactions
        for m, p in zip(['m1', 'm2', 'm3'], ['x2', 'x3', 'x1']):
            Protein.c(p) >> Protein.c(p) + mRNA.c(m)[lambda pro: f'{beta_m}/(1 + ({pro}/{k})^{n})']

        # Production reactions
        for m, p in zip(['m1', 'm2', 'm3'], ['x1', 'x2', 'x3']):
            mRNA.c(m) >> mRNA.c(m) + Protein.c(p)[beta_p]

        # We need the rate of degradation to be different from proteins and mRNA
        Mortal >> Zero[lambda r1: gamma_p if r1.is_a(Protein) else gamma_m]
        # This is the leaky mRNA expression, it needs to be low
        Zero >> Creator[leaky]

        MySim = Simulation(mRNA | Protein)
        return MySim.compile()

    assert compare_model(oscillator(), 'test_tools/model_7.txt')


# Model to test well defined orthogonal spaces
def orthogonal_spaces():
    try:
        A, B = BaseSpecies(2)
        A.a, A.b
        C = New(B)
        C.a, C.b
        MySim = Simulation(A | C)
        MySim.level = 0
        MySim.compile()
        return False
    except SystemExit:
        return True


def test_orthogonal():
    assert orthogonal_spaces()


# Model to test dimensional inconsistency
def dimensional_inconsistency():
    try:
        A, B, C = BaseSpecies(3)
        A(1 * u.mol / u.meter ** 3) + B(1 * u.mol / u.meter ** 2) >> C[1]
        MySim = Simulation(A | B | C)
        MySim.level = 0
        MySim.compile()
        return False
    except SystemExit:
        print('Dimensional inconsistency model Ok')
        return True


def test_dimensional_inconsistency():
    assert dimensional_inconsistency()


def average_value():
    try:
        E = BaseSpecies(1)
        Zero >> E[12]
        E >> Zero[25]

        MySim = Simulation(E)
        MySim.save_data = False
        MySim.run()
    except:
        pass


def hybrid_sim():
    A, B = BaseSpecies(2)
    A >> 2 * A[1]

    A(1)
    S1 = Simulation(A | B)
    S1.save_data = False
    S1.plot_data = False
    S1.duration = 3

    A.reset_reactions()
    A + B >> Zero[0.01]

    B(50)
    S2 = Simulation(A | B)
    S2.method = 'stochastic'
    S2.duration = (A <= 0) | (B <= 0)

    Sim = S1 + S2
    return compare_model(Sim.compile(), 'test_tools/model_8.txt')


def test_hybrid_sim():
    assert hybrid_sim()


def concatenated_simulation():
    A, B, C = BaseSpecies(3)
    A >> Zero[1]

    A(50)
    S1 = Simulation(A)
    S1.plot_data = False
    S1.duration = 5

    B >> Zero[1]

    B(50)
    S2 = Simulation(B)
    S2.duration = 5

    C >> Zero[1]

    C(50)
    S3 = Simulation(C)
    S3.duration = 5

    S = S1 + S2 + S3
    S.run()
    return S.results[A][-1] < 1 and S.results[B][-1] < 1 and S.results[C][-1] < 1


def test_concatenated_simulation():
    assert concatenated_simulation()


def event_type_test():
    A, B, C, D, E, F = BaseSpecies(6)

    A + B >> Zero[1]

    A(50), B(50), C(0)
    S = Simulation(A | B | C | D | E | F)
    S.plot_data = False

    with S.event_time(0) as _:
        F(1)

    with S.event_condition() as _:
        if (A <= 1) & (B <= 1):
            C(1)
        if A <= 1:
            D(1)
        if B <= 1:
            E(1)

    S.duration = 5
    S.run()
    return S.results[C][-1] == 1 and S.results[D][-1] == 1 and S.results[E][-1] == 1 and S.results[F][-1] == 1


def test_event_type_test():
    event_type_test()


def reacting_species_event():
    B = BaseSpecies(1)
    B.b1, B.b2
    A = New(B)

    B >> Zero[1]
    A.b2 >> Zero[0.5]
    A.a1 >> Zero[1]

    A.a1(100), A.b2(100), B.b1(100)
    S = Simulation(A | B)

    with S.event_condition() as _:
        if (A.a1 <= 10) & (B.b1 <= 10):
            A.a1(100)

    S.duration = 5
    return compare_model(S.compile(), 'test_tools/model_9.txt')


def test_reacting_species_event():
    assert reacting_species_event()


def unit_event_test():
    A = BaseSpecies(1)
    A >> Zero[1 / u.s]

    A(1 * u.mol)
    S = Simulation(A)
    with S.event_condition() as _:
        if A < 0.5 * u.mol:
            A(1 * u.mol)
    S.duration = 3
    return compare_model(S.compile(), 'test_tools/model_10.txt')


def test_unit_event_test():
    assert unit_event_test()


def reaction_deactivation():
    A, R = BaseSpecies(2)
    A + R >> 2 * A + R[1]

    A(1), R(1)
    S1 = Simulation(A | R)
    S1.duration = 1
    S1.plot_data = False

    R(0)
    S2 = Simulation(A | R)
    S2.duration = 1

    Sim = S1 + S2
    Sim.run()
    assert Sim.results[A][0] < Sim.results[A][-1] and Sim.results[R][0] == 1 and Sim.results[R][-1] == 0


def count_assignment():
    A = BaseSpecies(1)
    B = New(A)
    A.a1, A.a2

    B.b1 >> Zero[1]

    A.a1(100), A.a2(100)
    B.b1(100), B.b2(100)
    S = Simulation(A | B)
    S.plot_data = False
    S.duration = 5
    S.run()
    return compare_model(S.compile(), 'test_tools/model_11.txt') \
           and 150 > S.results[B][-1] > 100 and S.results[A][-1] == 200


def test_count_assignment():
    assert count_assignment()


def reaction_deactivation():
    A, R = BaseSpecies(2)
    A + R.r1 >> 2 * A + R.r1[1]

    A(1), R(1)
    S1 = Simulation(A | R)
    S1.plot_data = False
    S1.duration = 2

    R(0)
    S2 = Simulation(A | R)
    S2.duration = 2
    Sim = S1 + S2
    Sim.run()
    return Sim.results[R][-1] == 0 and Sim.results[R][0] == 1 and Sim.results[A][-2] == Sim.results[A][-1]


def test_reaction_deactivation():
    assert reaction_deactivation()


def complex_cell_model():
    Resource, Phage, Infectable = BaseSpecies(3)
    Cell = New(Infectable)
    Cell.t1, Cell.t2, Cell.t3

    def reproduction_rate(r):
        if r.t1:
            return 2
        elif r.t2:
            return 1
        elif r.t3:
            return 0.5

    Infectable.not_infected + Phage >> Infectable.infected[1]
    Cell + Resource >> 2 * Cell[lambda r1: reproduction_rate(r1)]
    Zero >> Resource[1]
    Resource >> Zero[1]
    Cell >> Zero[0.1]

    Cell.t1(1), Cell.t2(1), Cell.t3(1)
    S1 = Simulation(Cell | Resource | Phage)
    S1.duration = 30

    Cell.reset_reactions()
    Phage(1000)
    S2 = Simulation(Cell | Phage)
    S2.duration = 10
    S = S1 + S2
    S.plot_data = False
    S.run()

    for i, c in enumerate(S.results['Cell.t1.not_infected']):
        if c > 0:
            continue
        else:
            change_index = i
            break

    boll_1 = round(S.results[Cell.t1.not_infected][change_index - 1], 2) \
             == round(S.results[Cell.t1.infected][-1], 2)
    boll_2 = round(S.results[Cell.t1.not_infected][change_index - 1], 2) \
             == round(S.results[Cell.t1.infected][-1], 2)
    boll_3 = round(S.results[Cell.t1.not_infected][change_index - 1], 2) \
             == round(S.results[Cell.t1.infected][-1], 2)
    return boll_1 and boll_2 and boll_3


def test_complex_cell_model():
    assert complex_cell_model()


def zero_rate_reactions():
    A, B = BaseSpecies(2)

    A.a1, A.a2, B.b1, B.b2
    Combination = A * B

    Combination >> Zero[lambda r1: 0 if r1.b2 else 1]
    S = Simulation(Combination)
    return compare_model(S.compile(), 'test_tools/model_12.txt')


def test_zero_rate_reaction():
    zero_rate_reactions()


def double_rate():
    A, B = BaseSpecies(2)
    A.a1, A.a2, B.b1, B.b2

    def rate(r1, ball):
        factor1 = 0.5 if r1.a1 else 1
        factor2 = 0.5 if ball.b1 else 1
        return factor1 * factor2

    A + B >> Zero[rate]
    S = Simulation(A | B)
    return compare_model(S.compile(), 'test_tools/model_13.txt')


def test_double_rate():
    assert double_rate()


def single_rate():
    A, B = BaseSpecies(2)
    A.a1, A.a2, B.b1, B.b2

    def rate(r1):
        factor1 = 0.5 if r1.a1 else 1
        return factor1

    A + B >> Zero[rate]
    S = Simulation(A | B)
    return compare_model(S.compile(), 'test_tools/model_14.txt')


def test_single_rate():
    assert single_rate()


def triple_rate():
    A, B = BaseSpecies(2)
    A.a1, A.a2, B.b1, B.b2

    def rate(r1, r2, r3):
        factor1 = 0.5 if r1.a1 else 1
        factor2 = 0.5 if r2.b1 else 1
        factor3 = 0.5 if r3.c1 else 1
        return factor1*factor2*factor3

    A + B >> Zero[rate]
    S = Simulation(A | B)
    return compare_model(S.compile(), 'test_tools/model_13.txt')


def test_triple_rate():
    assert triple_rate()
