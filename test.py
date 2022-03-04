from mobspy import *
import os


if __name__ == '__main__':

    A, B = BaseSpecies(2)
    A(10)
    Zero >> B [lambda: f'{A}']
    Simulation(A | B).compile()

    exit()

    A, B, C, D = BaseSpecies(4)
    Test = New(A, 1)

    Test.b(30) >> B [lambda: f'3*1/{Test.c}']
    Test.c(30) >> C [lambda: f'3*1/{Test.c}']
    Test.d(30) >> D [lambda: f'3*1/{Test.d}']
    MySim = Simulation(Test | B | C | D)
    MySim.save_data = False
    MySim.run()
    exit()


