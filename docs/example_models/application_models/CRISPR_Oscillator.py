import sys, os
from mobspy import *

"""
    This is the CRISP oscillator from the MobsPy publication
    We start by declaring 3 base species and assigning their characteristics
    We need introduce characteristics here for the queries in other dimensions later
"""
Promoter, dCas, CasBinding = BaseSpecies()
Promoter.active, Promoter.inactive, CasBinding.no_cas, CasBinding.cas

"""
    The New function creates a New species from Promoter and CasBinding
    The created species inherit the characteristics from Promoter and CasBinding
    Thus allowing us to reference the characteristics from the old species and also
    add a new set of independent characteristics to them
"""
DNAPro = New(Promoter)
gRNA = New(CasBinding)

"""
    The reason why the New constructor was used is to add the list of different gRNAs and Promoters
    If they are added to the DNAPro or gRNA species, those states become independent from 
    active, inactive and no_cas, cas respectively 
"""
G = ListSpecies(3, gRNA)
P = ListSpecies(3, Promoter)

rev_rt = (1.8e-3/(u.nanomolar * u.second), 2.3e-2 / u.minute)
Rev[gRNA.no_cas + dCas >> gRNA.cas][rev_rt]
gRNA.no_cas >> Zero[0.0069 / u.second]

Promoter.active >> 2 * Promoter.active[2.3e-2 / u.minute]
Promoter >> Zero[2.3e-2 / u.minute]

"""
    Simple loop through the gRNAs and Promoters for assigning the activation reaction
"""
for Prom, Grna in zip(P, G):
    act_rt = lambda dna: 5 / u.minute if dna.active else 0
    Prom >> Grna.no_cas + Prom[act_rt]

"""
    Simple loop for the repression reaction
    Here we use the characteristics of gRNA_rep_list as the different type of gRNAs
"""
gRNA_rep_List = [G[-1], G[0], G[1]]
for Prom, Grna in zip(P, gRNA_rep_List):
    dna_rt1 = 1.2e-2 * u.liter / (u.nanomoles * u.second)
    dna_rt2 = 2.3e-2 / u.minute
    Prom.active + Grna.cas >> Prom.inactive [dna_rt1]
    Prom.inactive >> 2 * Prom.active + Grna.cas [dna_rt2]

# Rev defines a reversible reaction in both senses
Rev[Zero >> dCas][1 / u.minute, 2.3e-2 / u.minute]

for i in range(3):
    P[i].active(1)
G[0].no_cas(0), G[1].no_cas(40*u.nanomolar), G[2].no_cas(1*u.nanomolar)
dCas(43*u.nanomolar)

MySim = Simulation(G | P | dCas)
MySim.volume = 1*u.femtoliter

# 40 nanomolar, 1
# dCas 43nm
# MySim.plot_data = False
MySim.duration = 650*u.hours
MySim.step_size = 100
MySim.unit_x = 'hours'
MySim.run()

