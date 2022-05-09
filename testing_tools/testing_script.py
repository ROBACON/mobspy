# This script tests the model creation and compilation
# It also test the calculation capabilities
# It uses some simple model and assertions
#

from mobspy import *
import sys

# TODO Plot has random order for species names

if __name__ == '__main__':

    # Compare results with expected file
    def compare_model(comp_results, file_name):

        with open(file_name, 'r') as file:
            for r, line in zip(comp_results.split('\n'), file.readlines()):
                line = line.replace('\n', '')
                assert r == line, f'{r} and {line} do not match'

        print(f'Model {file_name} passed the test', file=sys.stderr)


    # Model to test the basics
    A, B, C = BaseSpecies(3)
    A + B >> C[1]
    MySim = Simulation(A | B | C)
    MySim.level = 0
    results = MySim.compile()
    compare_model(results, 'model_1.txt')

    # Model to test basic inheritance
    Carnivore, Herbivore = BaseSpecies(2)
    Cat, Dog = New(Carnivore, 2)
    Carnivore + Herbivore(1 * u.mol) >> Carnivore[1]
    Cat(1 * u.mol), Dog(1 * u.mol)
    MySim = Simulation(Cat | Dog | Herbivore)
    MySim.level = 0
    MySim.volume = 1 * u.meter ** 2
    results = MySim.compile()
    compare_model(results, 'model_2.txt')

    # Model to test species multiplication
    MGMT, Blue_Oyster_Cult, The_Smiths = BaseSpecies(3)
    MGMT.eletric_fell, MGMT.little_dark_age, MGMT.kids
    Blue_Oyster_Cult.burning_for_you >> Blue_Oyster_Cult.reaper[1]
    The_Smiths.stop_me >> The_Smiths.charming_man[1]
    Music = MGMT * Blue_Oyster_Cult * The_Smiths
    MySim = Simulation(Music)
    MySim.level = 0
    results = MySim.compile()
    compare_model(results, 'model_3.txt')

    # Model to test inheritance queries
    # All bacterias are infected by any virus here
    Bacteria, Virus = BaseSpecies(2)
    B1, B2 = New(Bacteria, 2)
    V1, V2 = New(Virus, 2)
    Bacteria.not_infected + Virus >> Bacteria.infected[1]
    MySim = Simulation(B1 | B2 | V1 | V2)
    MySim.level = 0
    results = MySim.compile()
    compare_model(results, 'model_4.txt')

    # Model to test round-robin and stoichiometry
    A = BaseSpecies(1)
    B, C = New(A, 2)
    A >> 2 * A[1]
    2 * A >> 3 * A[1]
    MySim = Simulation(B | C)
    MySim.level = 0
    MySim.compile()
    compare_model(results, 'model_5.txt')

    # Model to test well defined orthogonal spaces
    try:
        A, B = BaseSpecies(2)
        A.a, A.b
        C = New(B)
        C.a, C.b
        MySim = Simulation(A | C)
        MySim.compile()
    except SystemExit:
        print('Different species with shared-characteristics error Ok')

    # Model to test dimensional inconsistency
    try:
        A, B = BaseSpecies(2)
        A(1 * u.mol / u.meter ** 3) + B(1 * u.mol / u.meter ** 2) >> C[1]
        MySim = Simulation(A | B | C)
        MySim.level = 0
        MySim.compile()
    except SystemExit:
        print('Dimensional inconsistency model Ok')






