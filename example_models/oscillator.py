import sys, os
import matplotlib.pyplot as plt
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *

# TODO look at repression and promoter reactions

Promoter, Chemical = BaseSpecies(2)
Promoter.inactive, Promoter.active
DNAPromoter = New(Promoter)
DNAPromoter.PLcl(100), DNAPromoter.PLacl(100), DNAPromoter.PTetR(100)
Chemical.TetR(1)

list_of_chemicals = ['TetR', 'Lcl', 'Lacl']
list_of_promoters = ['PLcl', 'PLacl', 'PTetR']
for che, pro in zip(list_of_chemicals, list_of_promoters):
    Rev[ DNAPromoter.inactive.v(pro) + Chemical.v(che) >> DNAPromoter.active.v(pro)][1, 1]
Chemical >> Zero [1]

repressed = ['TetR', 'Lcl', 'Lacl']
repressors = ['Lacl' ,'TetR', 'Lcl']
hill = lambda che: f'10/(1 + ({che})^3)'
for rpsor, rpsed in zip(repressed, repressors):
    Chemical.v(rpsor) >> Chemical.v(rpsed) + Chemical.v(rpsor) [hill]

MySim = Simulation(DNAPromoter | Chemical)
MySim.duration = 300
MySim.plot_data = False
MySim.save_data = False
MySim.run()
MySim.plot_deterministic(Chemical.TetR, Chemical.Lcl, Chemical.Lacl)



