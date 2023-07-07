from mobspy import *
import numpy as np

Activatable, Age, Infected, Mortal, Signal, Phage = BaseSpecies()
Signal, Phage = New(Mortal)
Reproducer = New(Age)

Mortal >> Zero [0.1/u.h]

Signal + Activatable.inactive >> Activatable.active [20/u.h]
Activatable.active >> Activatable.inactive + Signal [0.01/u.h]

Age.young >> Age.old [1/u.h]

Reproducer >> 2*Reproducer.young [0.1/u.h]

infection_rate = lambda r1, r2: 2/u.h if r1.old else 1/u.h
Infected.not_infected + Phage >> Infected.infected [infection_rate]

Cell = Activatable*Reproducer*Infected*Mortal

Cell.active >> Cell.active + Signal [10/u.h]
Cell.active.infected >> Cell.active.infected + Phage [10/u.h]

x = [float(x) for x in np.linspace(0, 15, 50)]
initial_signal = ModelParameters(x)

counts = {Cell.not_infected: 50, Cell.infected: 50, Phage: 0, Signal: initial_signal}
model = set_counts(counts)
S = Simulation(model)
S.duration = 0.3*u.h
S.plot_data = False
S.unit_x = 'hours'
S.run()

S.plot(Cell.not_infected, Cell.infected)