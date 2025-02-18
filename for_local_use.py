from mobspy import *
from mobspy import modules

if __name__ == '__main__':

    L, R = BaseSpecies()

    kf, kr = ModelParameters(1e-3, 1e-3)
    L.sl_0 + R.sr_0 >> L.sl_1 + R.sr_1[kf, lambda r: kr * r]

    S = Simulation(L | R)
    print(S.compile())



