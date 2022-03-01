import sys, os
import matplotlib.pyplot as plt
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *

# TODO add bind function

Chemical, Promoter, Cell, Output = BaseSpecies(4)

Ara = Chemical*New
aTc = Chemical*New
Lasl = Chemical*New
Rhll = Chemical*New
Pbad = Promoter*New
Ptet = Promoter*New
Plas = Promoter*New
Prhll = Promoter*New
Promoter.inactive, Promoter.active


for comp in ['c1', 'c2']:
    Ara + Cell.c(comp) + Pbad.inactive >> Cell.c(comp) + Pbad.active [1]
for comp in ['c1', 'c3']:
    aTc + Cell.c(comp) + Ptet.inactive >> Cell.c(comp) + Ptet.active [1]
for comp in ['c2', 'c3']:
    Lasl + Cell.c(comp) + Plas.inactive >> Cell.c(comp) + Plas.active [1]
Rhll + Cell.c4 + Prhll.inactive >> Cell.c4 + Prhll.active [1]

Pbad.inactive + Ptet.inactive >> Lasl [1]

Pbad.inactive + Plas.inactive >> Rhll [0.01]
Plas.inactive + Ptet.inactive >> Rhll [0.01]

Prhll.active >> Output [0.01]

Cell.c1(1)
Cell.c2(1)
Cell.c3(1)
Cell.c4(1)
Pbad(10)
Ptet(10)
Plas(10)
Prhll(10)
MySim = Simulation(Cell | Ara | aTc | Lasl | Pbad | Ptet | Plas | Rhll | Prhll | Output)
MySim.save_data = False
MySim.plot_data = False
MySim.duration = 300
MySim.run()
MySim.plot_deterministic(Ara,aTc,Rhll,Output)

