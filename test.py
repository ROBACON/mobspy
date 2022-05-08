from mobspy import *

# TODO Plot has random order for species names

if __name__ == '__main__':

    A, B, C = BaseSpecies(3)
    A + B >> C[1]
    MySim = Simulation(A | B | C)
    results = MySim.compile()
    print('Bellow This:')
    print(results)
    print('Ends Here')


