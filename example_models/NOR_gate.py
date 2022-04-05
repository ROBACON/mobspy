import plotly.express as px
from mobspy import *


def NOR_GATE(A_conc, B_conc):
    A, B, Promoter, Protein = BaseSpecies(4)
    A + Promoter.inactive >> Promoter.active [0.5]
    B + Promoter.inactive >> Promoter.active[0.5]
    Promoter >> Promoter + Protein [lambda promoter: 1 if promoter.inactive else 0]
    Protein >> Zero [2]

    Promoter(100)
    A(A_conc), B(B_conc)
    MySim = Simulation(A | B | Promoter | Protein)
    MySim.duration = 10
    MySim.save_data = False
    MySim.plot_data = False
    MySim.run()
    return MySim.results['data']['Protein']['runs'][0][-1]


heatmap = []
for a in [0, 25, 50, 75, 100]:
    heatmap_line = []
    for b in [0, 25, 50, 75, 100]:
        output = NOR_GATE(a, b)
        heatmap_line.append(output)
    heatmap.append(heatmap_line)

for line in heatmap:
    print(line)

fig = px.imshow(heatmap, x=[0, 25, 50, 75, 100], y=[0, 25, 50, 75, 100], labels=dict(x='b', y='a'))
fig.show()