# This script tests the model creation and compilation
# It also test the calculation capabilities
# It uses some simple model and assertions

import pytest
from mobspy import *
import sys

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
