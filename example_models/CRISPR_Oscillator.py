#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:42:50 2022

@author: gayathriprakash
"""
import sys, os
import matplotlib.pyplot as plt

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

Rev[ gRNA.no_cas + dCas >> gRNA.cas ][0.0018*6.022e-1, 0.000383]
gRNA.no_cas >> Zero [0.0069]

Promoter.active >> 2*Promoter.active [0.000383]
Promoter >> Zero [0.000383]
for prom, grna in zip(Promoters_list,gRNAs_list):
    DNAPromoter.c(prom) >> gRNA.no_cas.c(grna) + DNAPromoter.c(prom)[lambda dna: 5/60 if dna.active else 0]
    
gRNA_rep_List = ['g3','g1','g2']
for prom, grna in zip(Promoters_list,gRNA_rep_List):
    DNAPromoter.active.c(prom) + gRNA.cas.c(grna) >> DNAPromoter.inactive.c(prom) [0.012*6.022e-1]
    DNAPromoter.inactive.c(prom) >> 2*DNAPromoter.active.c(prom) + gRNA.cas.c(grna) [0.000383]
    
Rev[ Zero >> dCas ] [1/60, 0.000383]

DNAPromoter.active.P1(1), DNAPromoter.active.P2(1), DNAPromoter.active.P3(1)
gRNA.no_cas.g1(0), gRNA.no_cas.g2(40/6.022e-1), gRNA.no_cas.g3(1/6.022e-1)
dCas(10)

MySim = Simulation( gRNA | DNAPromoter | dCas )
MySim.save_data = False
MySim.plot_data = False
MySim.duration = 8e5
MySim.output_event = True
MySim.step_size = 100

MySim.run()
MySim.plot_deterministic(gRNA.g1, gRNA.g2, gRNA.g3)

