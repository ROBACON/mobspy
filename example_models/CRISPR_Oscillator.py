import sys, os
from mobspy import *

Promoter, dCas, CasBinding = BaseSpecies(3)
Promoter.active, Promoter.inactive, CasBinding.no_cas, CasBinding.cas

DNAPro = New(Promoter)
gRNA = New(CasBinding)

gRNAs_list = ['g1', 'g2', 'g3']
Promoters_list = ['P1', 'P2', 'P3']

rev_rt = (1.8e-3/(u.nanomolar * u.second), 2.3e-2 / u.minute)
Rev[gRNA.no_cas + dCas >> gRNA.cas][rev_rt]
gRNA.no_cas >> Zero[0.0069 / u.second]

Promoter.active >> 2 * Promoter.active[2.3e-2 / u.minute]
Promoter >> Zero[2.3e-2 / u.minute]

for prom, grna in zip(Promoters_list, gRNAs_list):
    act_rt = lambda dna: 5 / u.minute if dna.active else 0
    DNAPro.c(prom) >> gRNA.no_cas.c(grna) + DNAPro.c(prom)[act_rt]

gRNA_rep_List = ['g3', 'g1', 'g2']
for prom, grna in zip(Promoters_list, gRNA_rep_List):
    dna_rt1 = 1.2e-2 * u.liter / (u.nanomoles * u.second)
    dna_rt2 = 2.3e-2 / u.minute
    DNAPro.active.c(prom) + gRNA.cas.c(grna) >> DNAPro.inactive.c(prom)[dna_rt1]
    DNAPro.inactive.c(prom) >> 2 * DNAPro.active.c(prom) + gRNA.cas.c(grna)[dna_rt2]

Rev[Zero >> dCas][1 / u.minute, 2.3e-2 / u.minute]

DNAPro.active.P1(1), DNAPro.active.P2(1), DNAPro.active.P3(1)
gRNA.no_cas.g1(0), gRNA.no_cas.g2(40*u.nanomolar), gRNA.no_cas.g3(1*u.nanomolar)
dCas(43*u.nanomolar)

MySim = Simulation(gRNA | DNAPro | dCas)
MySim.volume = 1*u.femtoliter

# 40 nanomolar, 1
# dCas 43nm
MySim.save_data = False
MySim.plot_data = False
MySim.duration = 2.5e6
MySim.output_event = True
MySim.step_size = 100

MySim.run()
MySim.plot_deterministic(gRNA.g1, gRNA.g2, gRNA.g3)

