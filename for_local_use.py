from mobspy import *
from mobspy import modules

if __name__ == '__main__':
    from mobspy import *


    # Replace parameters using units
    Color, Location = BaseSpecies()
    Color.red, Color.blue
    Location.here, Location.there
    Something = Color*Location

    2*Something >> 3*Something [lambda r1, r2: 1*u.decimeter**2/u.h if Location(r1) == Location(r2) else 0.5*u.decimeter**2/u.h]

    S = Simulation(Something)
    S.volume = 1*u.m**2
    print(S.compile())

    # print(A.get_characteristics())


