import sys, os
import matplotlib.pyplot as plt

# What does this part below do?
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *

"""
    This is an example of a full non-compacted model to help people understand MobsPy
    Compare it with CRISP_Oscillator
"""

# Adding Cas part to gRNA
Promoter, dCas, gRNA = BaseSpecies(3)


Promoter.active, Promoter.inactive
gRNA.no_cas, gRNA.cas
Zero >> dCas [1/60]
dCas >> Zero [0.000383]


g1, g2,g3 = New(gRNA,3)
P1,P2,P3 = New(Promoter,3)


g1.no_cas >> Zero[0.0069]
g2.no_cas >> Zero[0.0069]
g3.no_cas >> Zero[0.0069]

g2.cas >> g2.no_cas + dCas [0.000383]
g1.cas >> g1.no_cas + dCas [0.000383]
g3.cas >> g3.no_cas + dCas [0.000383]

g2.no_cas + dCas >> g2.cas [6.022e-1*0.0018]
g1.no_cas + dCas >> g1.cas [6.022e-1*0.0018]
g3.no_cas + dCas >> g3.cas [6.022e-1*0.0018]

P1.active >> g1.no_cas + P1.active [5/60]
P2.active >> g2.no_cas + P2.active [5/60]
P3.active >> g3.no_cas + P3.active [5/60]

P1.active + g3.cas >> P1.inactive [0.012*6.022e-1]
P2.active + g1.cas >> P2.inactive [0.012*6.022e-1]
P3.active + g2.cas >> P3.inactive [0.012*6.022e-1]


P1.inactive >> P1.active + P1.active + g3.cas [0.000383]
P2.inactive >> P2.active + P2.active + g1.cas [0.000383]
P3.inactive >> P3.active + P3.active + g2.cas [0.000383]


# DNA replication
P1.active >> P1.active + P1.active [0.000383]
P2.active >> P2.active + P2.active [0.000383]
P3.active >> P3.active + P3.active [0.000383]

# Promoter degradation
P1 >> Zero [0.000383]
P2 >> Zero [0.000383]
P3 >> Zero [0.000383]


P1.active (1)
P2.active (1)
P3.active (1)
g2.no_cas(40/6.022e-1)
g1.no_cas(0)
g3.no_cas(1/6.022e-1)
dCas(10)


MySim = Simulation(g1 | g2 | g3 | P1 | P2 | P3 | dCas )

MySim.save_data = False
MySim.plot_data = False
MySim.duration = 2000000


MySim.run()
#MySim.plot.logscale = ['Y']
MySim.plot_deterministic(g1,g2,g3)