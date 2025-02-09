from mobspy import *
import seaborn
import matplotlib.pyplot as plt
import numpy as np


max_plas, max_pbad, max_ptet  = (1, 1, 1)

Mortal, Location = BaseSpecies()
Location.c1, Location.c2, Location.c3, Location.c4
PBad, PTet, Pcl, PLas, Movable = New(Location)
Ara, aTc, Cl, YFP, AHL, OC12 = New(Mortal*Movable)
n = 50
# n = 10 - for testing
# Parametric sweeps for XOR gate analysis
ara_entrace_rate, atc_entrace_rate = ModelParameters([float(x) for x in np.linspace(0, 35, n)],
                                                     [float(x) for x in np.linspace(0, 2000, n)])

# Death and movement reactions
Mortal >> Zero[1]
for x, y in zip(['c1', 'c1', 'c2', 'c3'], ['c2', 'c3', 'c4', 'c4']):
    Movable.c(x) >> Movable.c(y) [0.1]

# Represents the entrance of aTc and Ara in the system
def diffusion_in_cell(Molecule, rate, locations):
    for l in locations:
        Zero >> Molecule.c(l) [rate]
diffusion_in_cell(Ara, ara_entrace_rate, ['c1', 'c2'])
diffusion_in_cell(aTc, atc_entrace_rate, ['c1', 'c3'])

# Promoter activation function
def promoter_activation(P, Ligand, Protein, tf_max, n,
                        K_d, protein_production_rate, locations):
    tf_linked = lambda lig: tf_max * lig ** n / (lig ** n + K_d ** n)
    tf_free = lambda lig: tf_max - tf_max * lig ** n / (lig ** n + K_d ** n)
    pr = lambda r1, r2: protein_production_rate(r1, tf_linked(r2), tf_free(r2))
    for l in locations:
        with Location.c(l):
            P + Ligand >> P + Ligand + Protein [pr]

# Inverter for each nor gate
def inverter_wire(P, R, Signal, locations):
    rate_f = lambda r1, r2: 181*r1*350/(1 + 350 + 15*r2 + 50*r2 + 15*50*0.18*r2**2)
    for l in locations:
        with Location.c(l):
            P + R >> P + R + Signal[rate_f]

# Custom buffer for clear visibility
def buffer(L, Signal, l, n, K):
    with Location.c(l):
        L >> L + Signal [lambda r: 30*r**n/(r**n + K**n)]

# Each promoter expression is written here and assign to the promoter function
pbad_p_rate = lambda r1, r2, r3: 765*r1*(0.009 + 37.5*r2)/(1 + 0.009 + 37.5*r2 + 3.4*r3)
promoter_activation(PBad, Ara, Cl, max_pbad, 2.8, 90, pbad_p_rate, ['c1', 'c2'])
ptet_p_rate = lambda r1, r2, r3: 300*r1*350/(1 + 350 + 2*160*r3 + 160**2*r3**2)
promoter_activation(PTet, aTc, Cl, max_ptet, 1.0, 250, ptet_p_rate, ['c1', 'c3'])
plas_p_rate = lambda r1, r2, r3: 69*r1*(0.002 + 100*r2)/(1 + 0.002 + 100*r2)
promoter_activation(PLas, AHL, Cl, max_plas, 1.4, 0.2, plas_p_rate, ['c2', 'c3'])
# Inverter wires in all the gates
inverter_wire(Pcl, Cl, AHL, ['c1']), inverter_wire(Pcl, Cl, OC12, ['c2', 'c3'])
# Buffer for YFP
buffer(OC12, YFP, 'c4', 4, 0.04)

# Count assignment and Sim
model = set_counts({Ara: 0, aTc: 0, Cl: 0, YFP: 0, PBad.c1: 1, PBad.c2: 1, PTet.c1: 1, PTet.c3: 1,
                    Pcl.c1: 1, Pcl.c2: 1, Pcl.c3: 1, PLas.c2: 1, PLas.c3: 1, AHL: 0, OC12: 0})
S = Simulation(model)
S.duration = 100
S.plot_data = False
S.run()

matrix = []
line = []
for i, r in enumerate(S.results[YFP]):
    line.append(r[-1])
    if (i + 1) % n == 0:
        matrix.append(line)
        line = []

ax = seaborn.heatmap(matrix)

# Invert the y-axis
ax.invert_yaxis()

plt.title('XOR Gate Heatmap', fontsize=18)
plt.xlabel('Ara (a.u.)', fontsize=14)
plt.ylabel('aTc (a.u.)', fontsize=14)
plt.show()
