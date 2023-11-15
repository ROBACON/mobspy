# This script tests the model creation and compilation
# It also test the calculation capabilities
# It uses some simple model and assertions

# import pytest
from mobspy import *
import sys
import matplotlib.pyplot as plt


# Compare results with expected file
def compare_model(comp_results, file_name):
    with open(file_name, 'r') as file:
        for r, line in zip(comp_results.split('\n'), file.readlines()):
            line = line.replace('\n', '')

            if r != line:
                print('file: ' + line)
                print('test: ' + r)
                return False
    return True


# Model to test the basics
def test_model_1():
    A, B, C = BaseSpecies(3)
    A + B >> C[1]
    MySim = Simulation(A | B | C)
    MySim.level = -1
    results = MySim.compile()
    assert compare_model(results, 'test_tools/model_1.txt')


# Model to test basic inheritance
def test_model_2():
    Carnivore, Herbivore = BaseSpecies(2)
    Cat, Dog = New(Carnivore, 2)
    Carnivore + Herbivore(1 * u.mol) >> Carnivore[1]
    Cat(1 * u.mol), Dog(1 * u.mol)
    MySim = Simulation(Cat | Dog | Herbivore)
    MySim.level = -1
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
    MySim.level = -1
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
    MySim.level = -1
    results = MySim.compile()
    assert compare_model(results, 'test_tools/model_4.txt')


# Model to test round-robin and stoichiometry
def test_model_5():
    A = BaseSpecies(1)
    B, C = New(A, 2)
    A >> 2 * A[1]
    2 * A >> 3 * A[1]
    MySim = Simulation(B | C)
    MySim.level = -1
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
    MySim.level = -1
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
        MySim.level = -1
        return MySim.compile()

    assert compare_model(oscillator(), 'test_tools/model_7.txt')


# Model to test well defined orthogonal spaces
def test_orthogonal_spaces():
    try:
        A, B = BaseSpecies(2)
        A.a, A.b
        C = New(B)
        C.a, C.b
        MySim = Simulation(A | C)
        MySim.level = -1
        MySim.compile()
        assert False
    except SystemExit:
        assert True


# Model to test dimensional inconsistency
def test_dimensional_inconsistency():
    try:
        A, B, C = BaseSpecies(3)
        A(1 * u.mol / u.meter ** 3) + B(1 * u.mol / u.meter ** 2) >> C[1]
        MySim = Simulation(A | B | C)
        MySim.level = -1
        MySim.compile()
        assert False
    except SystemExit:
        print('Dimensional inconsistency model Ok')
        assert True


def test_average_value():
    E = BaseSpecies(1)
    Zero >> E[12]
    E >> Zero[25]

    MySim = Simulation(E)
    MySim.save_data = False
    MySim.plot_data = False
    MySim.level = -1
    MySim.run()


def test_hybrid_sim():
    A, B = BaseSpecies(2)
    A >> 2 * A[1]

    A(1), B(50)
    S1 = Simulation(A)
    S1.save_data = False
    S1.plot_data = False
    S1.duration = 3
    S1.level = -1

    A.reset_reactions()
    A + B >> Zero[0.01]

    S2 = Simulation(A | B)
    S2.method = 'stochastic'
    S2.duration = (A <= 0) | (B <= 0)
    S2.level = -1

    Sim = S1 + S2
    Sim.run()

    assert compare_model(Sim.compile(), 'test_tools/model_8.txt')
    assert Sim.fres[A][-1] == 0 or Sim.fres[B][-1] == 0


def test_concatenated_simulation():
    A, B, C = BaseSpecies(3)
    A >> Zero[1]

    A(50)
    S1 = Simulation(A)
    S1.plot_data = False
    S1.duration = 5
    S1.level = -1

    B >> Zero[1]

    B(50)
    S2 = Simulation(B)
    S2.duration = 5
    S2.level = -1

    C >> Zero[1]

    C(50)
    S3 = Simulation(C)
    S3.duration = 5
    S3.level = -1

    S = S1 + S2 + S3
    S.run()
    assert S.fres[A][-1] < 1 and S.fres[B][-1] < 1 and S.fres[C][-1] < 1


def test_event_type():
    A, B, C, D, E, F = BaseSpecies(6)

    A + B >> Zero[1]

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
    assert compare_model(S.compile(), 'test_tools/model_15.txt')


def test_reacting_species_event():
    B = BaseSpecies(1)
    B.b1, B.b2
    A = New(B)

    B >> Zero[1]
    A.b2 >> Zero[0.5]
    A.a1 >> Zero[1]

    A.a1(100), A.b2(100), B.b1(100)
    S = Simulation(A | B)
    S.level = -1

    with S.event_condition((A.a1 <= 10) & (B.b1 <= 10)):
        A.a1(100)

    S.duration = 5
    assert compare_model(S.compile(), 'test_tools/model_9.txt')


def test_unit_event_test():
    A = BaseSpecies(1)
    A >> Zero[1 / u.s]

    A(1 * u.mol)
    S = Simulation(A)
    S.level = -1
    with S.event_condition(A < 0.5 * u.mol):
        A(1 * u.mol)
    S.duration = 3
    assert compare_model(S.compile(), 'test_tools/model_10.txt')


def test_reaction_deactivation():
    A, R = BaseSpecies(2)
    A + R >> 2 * A + R[1]

    # Adding this reaction for compatibility reasons with python versions lower than 3.10
    R >> Zero[1e-100]

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
    Sim.run()

    assert Sim.fres[A][0] < Sim.fres[A][-1] and Sim.fres[R][0] == 1 and Sim.fres[R][-1] == 0


def test_count_assignment():
    A = BaseSpecies(1)
    B = New(A)
    A.a1, A.a2

    B.b1 >> Zero[1]

    A.a1(100), A.a2(100)
    B.b1(100), B.b2(100)
    S = Simulation(A | B)
    S.level = -1
    S.plot_data = False
    S.duration = 5
    S.run()
    assert compare_model(S.compile(), 'test_tools/model_11.txt') \
           and 150 > S.fres[B][-1] > 100 and S.fres[A][-1] == 200


def test_complex_cell_model():
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
    S1.level = -1
    S1.duration = 30

    Cell.reset_reactions()
    Phage(1000)
    S2 = Simulation(Cell | Phage)
    S2.level = -1
    S2.duration = 10
    S = S1 + S2
    S.plot_data = False
    S.run()

    for i, c in enumerate(S.fres['Cell.t1.not_infected']):
        if c > 0:
            continue
        else:
            change_index = i
            break

    boll_1 = round(S.fres[Cell.t1.not_infected][change_index - 1], 2) \
             == round(S.fres[Cell.t1.infected][-1], 2)
    boll_2 = round(S.fres[Cell.t1.not_infected][change_index - 1], 2) \
             == round(S.fres[Cell.t1.infected][-1], 2)
    boll_3 = round(S.fres[Cell.t1.not_infected][change_index - 1], 2) \
             == round(S.fres[Cell.t1.infected][-1], 2)
    assert boll_1 and boll_2 and boll_3


def test_zero_rate_reactions():
    A, B = BaseSpecies(2)

    A.a1, A.a2, B.b1, B.b2
    Combination = A * B

    Combination >> Zero[lambda r1: 0 if r1.b2 else 1]
    S = Simulation(Combination)
    S.level = -1
    return compare_model(S.compile(), 'test_tools/model_12.txt')


def test_double_rate():
    A, B = BaseSpecies(2)
    A.a1, A.a2, B.b1, B.b2

    def rate(r1, ball):
        factor1 = 0.5 if r1.a1 else 1
        factor2 = 0.5 if ball.b1 else 1
        return factor1 * factor2

    A + B >> Zero[rate]
    S = Simulation(A | B)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_13.txt')


def test_single_rate():
    A, B = BaseSpecies(2)
    A.a1, A.a2, B.b1, B.b2

    def rate(r1):
        factor1 = 0.5 if r1.a1 else 1
        return factor1

    A + B >> Zero[rate]
    S = Simulation(A | B)
    S.level = - 1
    assert compare_model(S.compile(), 'test_tools/model_14.txt')


def test_triple_rate():
    A, B = BaseSpecies(2)
    A.a1, A.a2, B.b1, B.b2

    def rate(r1, r2, r3):
        factor1 = 0.5 if r1.a1 else 1
        factor2 = 0.5 if r2.b1 else 1
        factor3 = 0.5 if r3.c1 else 1
        return factor1*factor2*factor3

    A + B >> Zero[rate]
    S = Simulation(A | B)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_13.txt')


def test_stochastic_event_duration():

    A, B = BaseSpecies(2)
    A + B >> Zero[0.01]

    A(100), B(100)
    S1 = Simulation(A | B)
    S1.save_data = False
    S1.plot_data = False
    S1.method = 'stochastic'
    S1.duration = (A <= 0) | (B <= 0)
    S1.level = -1
    S1.run()
    R = S1.fres
    assert R[A][0] > 0 and R[B][0] > 0 and R[A][-1] == 0 and R[B][-1] == 0


def test_logic_operator_syntax():
    test_failed = False
    simlog.global_simlog_level = -1
    A, B = BaseSpecies(2)
    A.a1, A.a2, A.a3

    try:
        (10 >= A) >= 10
        test_failed = True
    except SystemExit:
        pass

    try:
        (10 >= A >= 10 >= A)
        test_failed = True
    except SystemExit:
        pass

    try:
        (10 >= A*A)
        test_failed = True
    except SystemExit:
        pass

    try:
        S1 = Simulation(A)
        S1.level = -1
        with S1.event_condition(B <= 10):
            A(100)
        S1.compile()
        test_failed = True
    except SystemExit:
        pass

    r1 = ((10 >= 2*A) & (A <= 10)) | (10 >= A)
    S = Simulation(A)
    S.level = -1
    with S.event_condition(r1):
        A(100)
    assert compare_model(S.compile(), 'test_tools/model_16.txt')

    if test_failed:
        assert False


def test_stack_position():

    Cell = BaseSpecies()
    Cell >> 2 * Cell[1]
    A, B, C = New(Cell)
    S = Simulation(A | B | C)
    S.level = -1
    compare_model(S.compile(), 'test_tools/model_17.txt')

    def hi():
        Cell = BaseSpecies()
        Cell >> 2 * Cell[1]
        A, B, C = New(Cell)

        S = Simulation(A | B | C)
        S.level = -1
        compare_model(S.compile(), 'test_tools/model_17.txt')
    hi()

    def hi_inside_hi():
        hi()
    hi_inside_hi()


def test_empty_arguments():
    A, B = BaseSpecies()

    A >> Zero[lambda: f'{A}*0.01']
    S = Simulation(A)
    S.duration = 5
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_18.txt')


def test_conditional_between_meta_species():
    Cu = BaseSpecies()
    Cu.c1, Cu.c2
    Azi, Byy = New(Cu)
    Azi.a1, Azi.a2, Byy.b1, Byy.b2

    Azi >> Zero[1]
    Byy >> Zero[0.1]

    Azi(200), Byy(50)
    S = Simulation(Azi | Byy)
    S.plot_data = False
    S.level = -1
    with S.event_condition(Azi.a1 <= Byy.b1):
        Azi(200)
    S.duration = 10
    S.method = 'stochastic'
    assert compare_model(S.compile(), 'test_tools/model_19.txt')
    S.run()
    for i in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
        assert S.fres[Azi][i] > S.fres[Byy][i]


def test_conditional_between_meta_species_2():
    A, B = BaseSpecies()

    A >> Zero[1]
    B >> Zero[0.1]

    r1 = ((A < B) & (A < B) | (A < B))

    A(200), B(50)
    S = Simulation(A | B)
    with S.event_condition(r1):
        A(200)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_20.txt')


def test_event_reaction_not_allowed():

    try:
        A = BaseSpecies()
        A >> Zero[1]

        S = Simulation(A)

        with S.event_time(0):
            Zero >> A[1]
        assert False
    except SystemExit:
        assert True


def all_test():
    A, B = BaseSpecies()
    A.a1, A.a2
    B.b1, B.b2

    C = A * B

    All[C](100)
    C >> All[C][1]

    S = Simulation(C)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_21.txt')


def all_test_2():
    B = BaseSpecies()
    B.b1, B.b2
    C, D = New(B)
    C.c1, C.c2, D.d1, D.d2

    Zero >> All[B.b1][1]

    S = Simulation(C | D)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_22.txt')


def test_error_mult():
    try:
        D = BaseSpecies(1)
        A, B, C = D * BaseSpecies(3)
        simlog.global_simlog_level = -1
        assert False
    except SystemExit:
        assert True


def test_set_counts():
    A, C = BaseSpecies()
    A.a1, A.a2
    B = New(A)
    B.b1, B.b2

    model = set_counts({All['B.a1']: 100, C:200*u.mols, 'A.a1': 100, A.a2: 50})
    S = Simulation(model)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_23.txt')


def test_bool_error():
    B = BaseSpecies()

    B >> Zero[1]
    B(100)

    S = Simulation(B)

    simlog.global_simlog_level = -1
    try:
        with S.event_condition(B == 0):
            B(100)
        assert False
    except SystemExit:
        pass

    try:
        S.duration = True
        assert False
    except SystemExit:
        pass
    assert True


def test_event_all():
    Acka = BaseSpecies()
    Acka.a1, Acka.a2

    Baka = New(Acka)
    Baka.b1, Baka.b2

    Baka >> Zero[1]

    S = Simulation(Baka)
    S.level = -1
    S.plot_data = False
    with S.event_time(0):
        set_counts({All[Baka]: 10, Baka.b1.a2: 20})
        Baka.a1.b1(30)
    assert compare_model(S.compile(), 'test_tools/model_24.txt')
    S.duration = 5
    S.step_size = 1
    S.run()
    assert S.fres[Baka.a1.b1][0] == 30
    assert S.fres[Baka.b1.a2][0] == 20
    assert S.fres[Baka.a1.b2][0] == 10
    assert S.fres['Baka.a1.b1'][0] == 30
    assert S.fres['Baka.b1.a2'][0] == 20
    assert S.fres['Baka.a1.b2'][0] == 10
    for key in S.fres:
        if key == 'Time':
            continue
        assert S.fres[key][-1] < 1


def test_one_value_concatenation_sim():
    A, B = BaseSpecies()

    B(200)
    S2 = Simulation(A | B)
    S2.plot_data = False
    S2.level = -1
    S2.duration = 5
    S2.step_size = 1
    S2.duration = (A <= 0) | (B <= 0)
    S2.run()
    assert len(S2.fres[A]) == 1


def test_crash_after_modification():

    try:
        A = BaseSpecies()
        S1 = Simulation(A)
        A = BaseSpecies()
        A.a1, A.a2

        S2 = Simulation(A)
        S = S1 + S2
        S.level = -1
        S.run()
        assert False
    except SystemExit:
        assert True


def test_unit_bi_dimension():
    A = BaseSpecies()
    A(5 / u.m ** 2)
    S = Simulation(A)
    S.volume = 2 * u.m ** 2
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_25.txt')


def test_bi_dimensional_rates():
    Ball, Child, Bacteria = BaseSpecies(3)

    Ball(10 / u.meter ** 2)
    Child(1 / u.meter ** 2)
    Bacteria(1 * u.mol)

    Bacteria >> Zero[1 * u.mol / u.second]
    Ball + Child + Child >> Ball + Child[1e-3 * (u.meter ** 4) / u.hour]
    Ball + Child >> Ball[1e-3 * (u.meter ** 2) / u.hour]

    S = Simulation(Ball | Child | Bacteria)
    S.volume = 2 * u.m ** 2
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_26.txt')


def test_dimension_in_function_only():
    A = BaseSpecies()

    A + A >> 3 * A[lambda: 1 * u.milliliter / u.second]

    A(1)
    S = Simulation(A)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_27.txt')


def test_multiple_simulation_counts():

    Age, Color, Size = BaseSpecies()

    Age.young, Age.old,
    Color.blue, Color.red,
    Size.small, Size.big

    Tree = Age * Color * Size

    Tree(100), Tree.red.big(150), All[Tree](10), All[Tree.big](100)
    S = Simulation(Tree)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_28.txt')

    Tree.reset_quantities()
    model = set_counts({All[Tree]: 30, 'Tree.blue.old': 100})
    S = Simulation(model)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_29.txt')


def test_string_events_assignment():

    A = BaseSpecies()

    A.a1, A.a2, A.a3

    S = Simulation(A)
    S.level = -1
    with S.event_time(5):
        All[A](f'{A} + 1')

    with S.event_time(10):
        All[A.a1](f'{A} + 1')

    with S.event_time(15):
        A.a1(f'{A} + 1')

    compare_model(S.compile(), 'test_tools/model_30.txt')


def test_plotting():
    Color, Disease = BaseSpecies()

    Color.blue, Color.red, Color.yellow
    Disease.not_sick, Disease.sick

    Disease.not_sick >> Disease.sick[1]

    Tree = Color * Disease

    Tree.yellow(20), Tree.red(20), Tree.blue(20)

    S = Simulation(Tree)
    S.level = -1
    S.method = 'stochastic'
    S.plot_data = False
    S.repetitions = 3
    S.step_size = 0.25
    S.duration = 3
    S.run()
    S.plot_config.save_to = 'test_plot_images/stochastic_tree.png'
    S.plot_stochastic(Tree.not_sick, Tree.sick)

    S.plot_config.save_to = 'test_plot_images/deterministic_tree.png'
    S.plot(Tree.not_sick, Tree.sick)

    S.plot_config.save_to = 'test_plot_images/constant_tree.png'
    S.plot()

    assert os.path.exists('test_plot_images/stochastic_tree.png')
    assert os.path.exists('test_plot_images/deterministic_tree.png')
    assert os.path.exists('test_plot_images/constant_tree.png')


def test_volume_after_sim():

    A = BaseSpecies()

    Zero >> A [42*1/(u.s*u.milliliter)]
    A >> Zero [1]

    S = Simulation(A)
    S.plot_data = False
    S.level = -1
    S.volume = 1 * u.milliliter
    S.run()
    assert int(S.fres[A][-1]) == 42


def order_model_str(data_for_sbml):
    species_for_sbml = data_for_sbml['species_for_sbml']
    mappings_for_sbml = data_for_sbml['mappings']
    parameters_for_sbml = data_for_sbml['parameters_for_sbml']
    reactions_for_sbml = data_for_sbml['reactions_for_sbml']
    events_for_sbml = data_for_sbml['events_for_sbml']

    model_str = '\n'
    model_str += 'Species' + '\n'
    species_alpha = list(sorted(species_for_sbml.keys()))
    for spe in species_alpha:
        model_str += spe.replace('_dot_', '.') + ',' + str(species_for_sbml[spe]) + '\n'

    model_str += '\n'
    model_str += 'Mappings' + '\n'
    mappings_alpha = list(sorted(mappings_for_sbml.keys()))
    for map in mappings_alpha:
        model_str += map + ' :' + '\n'
        for element in sorted(mappings_for_sbml[map]):
            model_str += element + '\n'

    model_str += '\n'
    model_str += 'Parameters' + '\n'
    parameters_alpha = list(sorted(parameters_for_sbml.keys()))
    for par in parameters_alpha:
        model_str += par + ',' + str(parameters_for_sbml[par][0]) + '\n'

    model_str += '\n'
    model_str += 'Reactions' + '\n'
    remove_phantom_reactions = deepcopy(reactions_for_sbml)
    to_remove = []
    for reaction in remove_phantom_reactions:
        if 'phantom' in reaction:
            to_remove.append(reaction)
    for r in to_remove:
        remove_phantom_reactions.pop(r, None)
    reaction_alpha = [str(x[1]).replace('_dot_', '.') for x in
                      list(sorted(remove_phantom_reactions.items(), key=lambda x: str(x[1])))]

    for i, reac in enumerate(reaction_alpha):
        model_str += 'reaction_' + str(i) + ',' + reac + '\n'

    if events_for_sbml != {}:
        model_str += '\n'
        model_str += 'Events' + '\n'
        list_to_sort = [str(events_for_sbml[key]) for key in events_for_sbml]
        list_to_sort = sorted(list_to_sort)
        for i in range(len(list_to_sort)):
            model_str += ('event_' + str(i) + ',' + list_to_sort[i] + '\n').replace('_dot_', '.')

    return model_str


def test_parameters_with_sbml():

    A = BaseSpecies()
    A.a1, A.a2
    a, b, c, d, f, h = ModelParameters([1, 2], [1, 2], [1, 2], [1, 2], [1, 2], [1, 2])

    A >> 2 * A[lambda: f'5*(b + c)/10']

    All[A](a)
    S1 = Simulation(A)

    with S1.event_time(0):
        A.a2(d)

    with S1.event_time(2):
        A.a1('a + b')

    with S1.event_time(f):
        A.a1(d)

    S1.duration = 3

    B = BaseSpecies()

    B >> 2 * B[h]

    B(a)
    S2 = Simulation(A | B)
    S2.duration = 2

    S = S1 + S2
    S.plot_data = False
    S.level = -1
    S.run()

    model_str = ''
    for parameter_sweep in S1.sbml_data_list:
        for data_for_sbml in parameter_sweep:
            model_str += order_model_str(data_for_sbml)

    assert compare_model(model_str, 'test_tools/model_31.txt')


def test_shared_parameter_name():
    try:
        A = BaseSpecies()
        a = ModelParameters([1, 2])
        a.rename('A')

        A >> 2 * A[a]

        set_counts({'A': a})
        S = Simulation(A)
        S.level = -1
        S.plot_data = False
        S.run()
        assert False
    except:
        assert True


def test_set_counts_parameters():
    A = BaseSpecies()
    a = ModelParameters([1, 2])

    A >> 2 * A[a]

    set_counts({'A': a})
    S = Simulation(A)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_32.txt')


def test_repeated_parameters():

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
    except SystemExit:
        assert True


def initial_expression_test():

    A, B, Hey = BaseSpecies()
    D = New(A)

    A >> 2 * A[lambda r: 1 / u.h * (1 + 10 / r)]
    A + B >> Zero[lambda r1, r2: (1 * u.millimolar / u.h) * (1 + 10 * u.millimolar / r1 + 20 * u.millimolar / r2)]
    Hey >> Zero[lambda r: 1 / u.h * (20 * r + 30 * r + 40 * r)]
    D >> 2 * D[lambda r: 20 / u.h * r]

    S = Simulation(A | B | Hey | D)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_33.txt')


def test_wrong_dimension_error():

    try:
        A, B = BaseSpecies()

        A >> 2*A [lambda r: 1/u.h*(1 + 10/u.decimeter**3/r)]

        S = Simulation(A)
        S.level = -1
        S.compile()
        assert False
    except SystemExit:
        assert True

    try:
        A, B = BaseSpecies()

        A >> 2*A [lambda r: (1/(u.h*u.decimeter**3))*(1 + 10/r)]

        S = Simulation(A)
        S.level = -1
        S.compile()
        assert False
    except SystemExit:
        assert True


def test_more_than_used():

    A = BaseSpecies()

    Zero >> A[lambda r1: 20]

    S = Simulation(A)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_34.txt')


def zero_rate_test():

    A = BaseSpecies()

    Zero >> A[0]

    S = Simulation(A)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_35.txt')


def test_wrong_rate():

    try:
        Ara, aTc = BaseSpecies()

        Ara >> 2 * Ara[aTc]

        S = Simulation(aTc | Ara)
        S.compile()
        assert False
    except SystemExit:
        assert True


def test_conversion_outside():
    n_0 = 10
    mu_g = 0.2 / u.h

    Cell, Lysis, AHL, LuxI = BaseSpecies()

    Cell >> 2 * Cell[lambda cell: mu_g * cell * (n_0 - cell)]

    MySim = Simulation(Cell)
    MySim.level = -1
    assert compare_model(MySim.compile(), 'test_tools/model_36.txt')


def test_first_characteristic_in_reacting_species():
    A = BaseSpecies()
    A.something
    B = New(A)

    for a in [1, 2, 3]:
        Zero >> B.something.c('at_' + str(a))[1]

    B(1)
    S = Simulation(B)
    S.level = -1
    assert compare_model(S.compile(), 'test_tools/model_37.txt')


def test_model_reference():

    Mortal = BaseSpecies()
    A, B = New(Mortal)

    A(100), B(200)
    S1 = Simulation(A | B)
    assert str(S1.model) == '[\'A\', \'B\']'

    S1.duration = 0.5

    C = New(Mortal)
    C(50)
    S2 = Simulation(A | B | C)
    assert str(S1.model) == '[\'A\', \'B\']'
    assert str(S2.model) == '[\'A\', \'B\', \'C\']'


def test_sbml_generation():

    A = BaseSpecies()

    A >> Zero[1]

    A(100)
    S = Simulation(A)
    S.level = -1
    text = ''
    for sbml in S.generate_sbml():
        text += sbml
    assert compare_model(text, 'test_tools/model_38.txt')


def test_multi_sim_sbml():

    A = BaseSpecies()

    A >> Zero[1]

    A(100)
    S1 = Simulation(A)
    S2 = Simulation(A)
    S = S1 + S2
    S.level = -1

    text = ''
    for sbml in S.generate_sbml():
        text += sbml
    assert compare_model(text, 'test_tools/model_39.txt')


def test_inline_comment():

    A = BaseSpecies()

    A >> Zero[1]  # Test comment
    assert True


def test_with_statement_any_and_species_characteristics():

    Age, Color, Dense = BaseSpecies()
    Age.old, Age.young
    Color.red, Color.green
    Dense.dense, Dense.sparse
    Tree = Age * Color * Dense
    Grass = Age * Color * Dense

    with Age.old, Dense.sparse :
        with Color.red :
            Tree >> Grass [1]
        with Color.blue :
            Tree >> Grass [1]
            Tree(10)
        Tree(9)
        All[Grass](1)
    with Any.young.green :
        Tree + Grass >> Tree + Tree [2]

    S1 = Simulation(Tree | Grass)
    S1.level = -1

    assert compare_model(S1.compile(), 'test_tools/model_40.txt')

    Age, Color, Dense = BaseSpecies()
    Age.old, Age.young
    Color.red, Color.green
    Dense.dense, Dense.sparse
    Tree = Age * Color * Dense
    Grass = Age * Color * Dense

    with Age.old, Dense.sparse :
        with Any.red:
                Tree >> Grass [1]
        with Color.blue :
            Tree >> Grass [1]
            Tree(10)
        Tree(9)
        All[Grass](1)

    S2 = Simulation(Tree | Grass)
    S2.level = -1

    assert compare_model(S2.compile(), 'test_tools/model_41.txt')


def test_with_statement_on_any_and_event():

    A = BaseSpecies()
    A.a1, A.a2
    S = Simulation(A)
    S.level = -1

    with Any.a2, S.event_condition(A <= 0):
        A(100)
        
    assert compare_model(S.compile(), 'test_tools/model_42.txt')

    


# This is here because pytest is slow - but this script works fine with pytest. Just make sure that the
# python version in terminal is above 3.10
test_list = [test_model_1, test_model_2, test_model_3, test_model_4, test_model_5, test_model_6, test_model_7,
             test_orthogonal_spaces, test_average_value, test_hybrid_sim, test_concatenated_simulation,
             test_event_type, test_reacting_species_event, test_unit_event_test, test_reaction_deactivation,
             test_double_rate, test_single_rate, test_triple_rate, test_stochastic_event_duration,
             test_logic_operator_syntax, test_stack_position, test_empty_arguments,
             test_conditional_between_meta_species, test_conditional_between_meta_species_2,
             test_event_reaction_not_allowed, all_test, all_test_2, test_error_mult, test_set_counts,
             test_bool_error, test_event_all, test_one_value_concatenation_sim, test_crash_after_modification,
             test_unit_bi_dimension, test_bi_dimensional_rates, test_dimension_in_function_only,
             test_multiple_simulation_counts, test_string_events_assignment, test_plotting,
             test_volume_after_sim, test_parameters_with_sbml, test_shared_parameter_name,
             test_set_counts_parameters, test_repeated_parameters, initial_expression_test,
             test_wrong_dimension_error, test_more_than_used, zero_rate_test, test_wrong_rate,
             test_conversion_outside, test_first_characteristic_in_reacting_species, test_model_reference,
             test_sbml_generation, test_multi_sim_sbml, test_inline_comment,
             test_with_statement_any_and_species_characteristics, test_with_statement_on_any_and_event]

sub_test = test_list


def perform_tests():

    any_failed = False
    for test in sub_test:
        try:
            test()
            print(f'Test {test} passed')
        except:
            print('\033[91m' + f'Test {test} failed'  + '\033[0m', file=sys.stderr)
            any_failed = True

    if any_failed:
        assert False


perform_tests()
