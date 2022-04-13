import sys, os
import matplotlib.pyplot as plt
from mobspy import *

# TODO look at repression and promoter reactions

"""
    This is a more classical oscillator
    Here the repression is implemented using a hill function as a string
    It's the classical circle of three compounds repressing each other
"""

Promoter, Chemical = BaseSpecies(2)
Promoter.inactive, Promoter.active
DNAPromoter = New(Promoter)
DNAPromoter.PLcl(100), DNAPromoter.PLacl(100), DNAPromoter.PTetR(100)
Chemical.TetR(1)

list_of_chemicals = ['TetR', 'Lcl', 'Lacl']
list_of_promoters = ['PLcl', 'PLacl', 'PTetR']
for che, pro in zip(list_of_chemicals, list_of_promoters):
    Rev[ DNAPromoter.inactive.c(pro) + Chemical.c(che) >> DNAPromoter.active.c(pro)][1, 1]
Chemical >> Zero [1]

repressed = ['TetR', 'Lcl', 'Lacl']
repressors = ['Lacl' ,'TetR', 'Lcl']
hill = lambda che: f'10/(1 + ({che})^3)'
for rpsor, rpsed in zip(repressed, repressors):
    Chemical.c(rpsor) >> Chemical.c(rpsed) + Chemical.c(rpsor) [hill]

MySim = Simulation(DNAPromoter | Chemical)
MySim.duration = 300
MySim.plot_data = False
MySim.save_data = False
MySim.run()
MySim.plot_deterministic(Chemical.TetR, Chemical.Lcl, Chemical.Lacl)



