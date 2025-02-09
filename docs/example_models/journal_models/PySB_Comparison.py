from mobspy import *

L, R = BaseSpecies()

L_0, R_0, kf, kr = ModelParameters(100, 200, 1e-3, 1e-3)
Rev[ L.sl_0(L_0) + R.sr_0(R_0) >> L.sl_1 + R.sr_1 ]  [kf, kr]

S = Simulation(L | R)
S.run(duration=100, plot_data=False)
S.plot(L.sl_1)

