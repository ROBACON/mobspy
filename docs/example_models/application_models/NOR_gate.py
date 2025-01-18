import plotly.express as px
from mobspy import *


"""
    Here we have a NOR_GATE 
    There are two possible repressors for the Promoter A and B
    If any of them bind to the Promoter the protein can no longer be expressed
"""


def NOR_GATE(A_conc, B_conc):
    # Here we define the Protein to be produced, the Promoter that will act as the gate
    # A and B are the inputs any of them can inactivate the Promoter so they inherit from Repressor
    Repressor, Promoter, Protein = BaseSpecies()
    Repressor + Promoter.active >> Promoter.inactive [0.5]
    A, B = New(Repressor)
    Promoter >> Promoter + Protein [lambda promoter: 1 if promoter.active else 0]
    Protein >> Zero [2]

    Promoter(100)
    A(A_conc), B(B_conc)
    MySim = Simulation(A | B | Promoter | Protein)
    MySim.duration = 10
    MySim.level = 0
    MySim.save_data = False
    MySim.plot_data = False
    MySim.run()
    return MySim.fres['Protein'][-1]


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