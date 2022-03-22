import sys, os

# What does this part below do?
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *

Promoter, dCas, CasBinding = BaseSpecies(3)
Promoter.active, Promoter.inactive, CasBinding.no_cas, CasBinding.cas

DNAPromoter = New(Promoter)
gRNA = New(CasBinding)

gRNAs_list = ['g1', 'g2', 'g3']
Promoters_list = ['P1', 'P2', 'P3']

# Rev[ gRNA.no_cas + dCas >> gRNA.cas ] [0.0018*6.022e-1/u.second, 2.3e-2/u.second]
Rev[gRNA.no_cas + dCas >> gRNA.cas][1.8e-3 * (u.liter / (u.nanomoles * u.second)), 2.3e-2 / u.minute]
gRNA.no_cas >> Zero[0.0069 / u.second]

Promoter.active >> 2 * Promoter.active[2.3e-2 / u.minute]
Promoter >> Zero[2.3e-2 / u.minute]

for prom, grna in zip(Promoters_list, gRNAs_list):
    DNAPromoter.c(prom) >> gRNA.no_cas.c(grna) + DNAPromoter.c(prom)[lambda dna: 5 / u.minute if dna.active else 0]

gRNA_rep_List = ['g3', 'g1', 'g2']
for prom, grna in zip(Promoters_list, gRNA_rep_List):
    DNAPromoter.active.c(prom) + gRNA.cas.c(grna) >> DNAPromoter.inactive.c(prom)[
        1.2e-2 * u.liter / (u.nanomoles * u.second)]
    DNAPromoter.inactive.c(prom) >> 2 * DNAPromoter.active.c(prom) + gRNA.cas.c(grna)[2.3e-2 / u.minute]

Rev[Zero >> dCas][1 / u.minute, 2.3e-2 / u.minute]

# Fix this!!!!!

DNAPromoter.active.P1(1), DNAPromoter.active.P2(1), DNAPromoter.active.P3(1)
gRNA.no_cas.g1(0), gRNA.no_cas.g2(40*u.nanomolar), gRNA.no_cas.g3(1*u.nanomolar)
dCas(43*u.nanomolar)

# 40 nanomolar, 1
# dCas 43nm

MySim = Simulation(gRNA | DNAPromoter | dCas)
MySim.volume = 1*u.femtoliter
MySim.save_data = False
MySim.plot_data = False
MySim.duration = 2.5e6
MySim.output_event = True
MySim.step_size = 100

MySim.run()
MySim.plot_deterministic(gRNA.g1, gRNA.g2, gRNA.g3)

