import sys, os
import matplotlib.pyplot as plt
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *


Chemical, Promoter, Cell, Output = BaseSpecies(4)
Promoter.inactive, Promoter.active
DNAPromoter = Promoter*New

List_Chemicals = ["Ara",  "Ara",  "aTc",  "aTc",   "Lasl", "Lasl", "Rhll"]
List_Colonies  = ["c1",   "c2",   "c1",   "c3",    "c2",   "c3",   "c4"]
List_Promoters = ["Pbad", "Pbad", "Ptet", "Pter",  "Plas", "Plas", "Prhl"]
for chem, colonie, prom in zip(List_Chemicals, List_Colonies, List_Promoters):
    Chemical.v(chem) + Cell.v(colonie) + DNAPromoter.inactive.v(prom) >> Cell.v(colonie) + DNAPromoter.active.v(prom) [1]

List_P1       = ["Pbad", "Pbad", "Ptet"]
List_P2       = ["Ptet", "Plas", "Plas"]
List_Chemical = ["Lasl", "Rhll", "Rhll"]
for prom1, prom2, chem in zip(List_P1, List_P2, List_Chemical):
    DNAPromoter.inactive.v(prom1) + DNAPromoter.inactive.v(prom2) >> Chemical.v(chem) [lambda p1, p2: 1 if p1.Pbad and p2.Ptet else 0.01]

DNAPromoter.active.Prhl >> Output [0.01]

Cell.c1(1), Cell.c2(1),Cell.c3(1), Cell.c4(1), DNAPromoter.Pbad(10), DNAPromoter.Ptet(10), DNAPromoter.Plas(10), DNAPromoter.Prhl(10)
MySim = Simulation(Chemical | DNAPromoter | Cell | Output)
MySim.save_data = False
MySim.plot_data = False
MySim.duration = 300
MySim.run()
MySim.plot_deterministic(Output)


