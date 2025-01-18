from mobspy import *
import seaborn
import matplotlib.pyplot as plt

A, B, C, Pa, Pb = BaseSpecies()
C >> Zero [1]


def Promoter_Rule(Promoter, Ligand, Protein, K):
    Promoter + Ligand >> Promoter + Ligand + Protein [lambda p, l: (p/u.h)*l**4/(K**4 + l**4)]


Promoter_Rule(Pa, A, C, 5)
Promoter_Rule(Pb, B, C, 8)


x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
initial_a = ModelParameters(x)
initial_b = ModelParameters(x)

model = set_counts({A: initial_a, B: initial_b, C: 0, Pa: 1, Pb: 1})
S = Simulation(model)
S.plot_data = False
S.duration = 30*u.h
S.run()

# Heat map construction
line = []
matrix = []
for i in range(len(S.results)):
    line.append(S.results[C][i][-1])

    if (i + 1) % len(x) == 0:
        matrix.append(line)
        line = []

ax = seaborn.heatmap(matrix)
ax.set(xlabel='A', ylabel='B')
plt.show()



