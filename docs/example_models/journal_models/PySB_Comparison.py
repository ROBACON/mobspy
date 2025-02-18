from mobspy import *

L, R = BaseSpecies()

L_0, R_0, kf, kr = ModelParameters(100, 200, 1e-3, 1e-3)
L.sl_0 + R.sr_0 >> L.sl_1 + R.sr_1 [kf, lambda r: kr*r]

L(L_0), R(R_0)
S = Simulation(L | R)
print(S.compile())

