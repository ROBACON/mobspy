import sys, os
import matplotlib.pyplot as plt

from mobspy import *

"""
    Here we have the implementation of a toggle switch
    Since MobsPy does not yet have an event implementation, we call the model again with changes made to it
    In the future we will implement a way to handle events, but for now events can be handled in a similar fashion
"""

p1 = []
p2 = []
T = []


# Toggle Switch
def set_model_initial_conditions(simulation, state_dict):
    simulation.compile(verbose=False)
    for key in state_dict:
        for spe in simulation._species_for_sbml:
            value = state_dict[key]
            species_state = key.replace('.', '_dot_')
            if species_state == spe:
                simulation._species_for_sbml[species_state] = value


def event_conditions(simulation, event_dict):
    for key in event_dict:
        for spe in simulation._species_for_sbml:
            value = event_dict[key]
            species_state = key.replace('.', '_dot_')
            if species_state == spe:
                simulation._species_for_sbml[species_state] = value


def toggle_switch(g1_count=0, g2_count=0, state_dict={}, duration=500000):
    Promoter, Cas_binding, dCas9 = BaseSpecies(3)
    Cas_binding.cas, Cas_binding.no_cas
    Promoter.inactive, Promoter.active
    gRNA = New(Cas_binding, 1)
    DNAPromoter = New(Promoter, 1)
    Promoter_list = ['P1', 'P2']
    gRNA_list = ['g1', 'g2']

    # g_RNA pdn
    for grna, prom in zip(gRNA_list, Promoter_list):
        DNAPromoter.active.c(prom) >> gRNA.no_cas.c(grna) + DNAPromoter.active.c(prom)[0.0833]

    # Degradation rxns
    gRNA.no_cas >> Zero[0.0069]
    Rev[Zero >> dCas9][1 / 60, 0.000383]

    # Promoter Pdn and Degradation
    Promoter.active >> 2 * Promoter.active[0.000383]
    Promoter >> Zero[0.000383]

    # gRNA:dCas9 Complex Formation
    Rev[gRNA.no_cas + dCas9 >> gRNA.cas][0.0018 * 6.022e-1, 0.000383]

    # Complex formation at the promoter and Complex falloff during Repn
    gRNA_rep_list = ['g2', 'g1']
    for grna, prom in zip(gRNA_rep_list, Promoter_list):
        DNAPromoter.active.c(prom) + gRNA.cas.c(grna) >> DNAPromoter.inactive.c(prom)[0.012 * 6.022e-1]
        DNAPromoter.inactive.c(prom) >> 2 * DNAPromoter.active.c(prom) + gRNA.cas.c(grna)[0.000383]

    MySim = Simulation(gRNA | DNAPromoter | dCas9)
    MySim.save_data = False
    MySim.plot_data = False
    MySim.duration = duration
    if not state_dict:
        DNAPromoter.P1.active(1), DNAPromoter.P2.active(1)
        gRNA.no_cas.g1(g1_count), gRNA.no_cas.g2(g2_count)
        dCas9(2)
        MySim.run()
    else:
        if g1_count != 0:
            event_dict = {'gRNA.no_cas.g1': g1_count}
        if g2_count != 0:
            event_dict = {'gRNA.no_cas.g2': g2_count}
        set_model_initial_conditions(MySim, state_dict)
        event_conditions(MySim, event_dict)
        MySim.run(compile=False)

    state_dict = {}
    for key in MySim.results['data']:
        if key == 'Time':
            pass
        else:
            state_dict[key] = MySim.results['data'][key]['runs'][0][-1]

    global p1
    global p2
    p1 += MySim.results['data']['DNAPromoter.P1.active']['runs'][0]
    p2 += MySim.results['data']['DNAPromoter.P2.active']['runs'][0]
    return state_dict


n = 1e10
state = toggle_switch(g1_count=10 / 6.022e-1, g2_count=0)
state = toggle_switch(g1_count=0, g2_count=n*10 / 6.022e-1, state_dict=state)
state = toggle_switch(g1_count=n*10/6.022e-1, g2_count=0, state_dict=state)
state = toggle_switch(g1_count=0, g2_count=n*10/6.022e-1, state_dict=state)
plt.plot(p1)
plt.plot(p2)
plt.show()

