"""Exercise the installed wheel with Python's isolated import mode (-I)."""

import math

from mobspy import BaseSpecies, Model, Simulation, Zero, u

A = BaseSpecies(["A"])
A >> Zero @ (1 / u.minute)
A(100)
simulation = Simulation(Model(A))
simulation.duration = 1 * u.minute
simulation.volume = 0.2 * u.mL
simulation.compile(verbose=False)
simulation.compile(verbose=False)
simulation.run(plot_data=False, level=0)
if not math.isclose(simulation.fres["A"][-1], 500 / math.e, rel_tol=1e-5):
    raise RuntimeError("Installed wheel failed the unit-aware decay smoke test")
print("Installed wheel smoke test passed")
