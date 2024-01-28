import plotly.express as px
from mobspy import *

"""
    In this model we have a simple AND_GATE
    Here we have the species A and B
    Once they encounter they form a compound AB that is able to activate a promoter
    With the promoter active a protein is produced
    The protein decays to zero (dies) at a certain rate
"""


def AND_GATE(A_conc, B_conc):
    # Here we define the Protein to be produced, the Promoter that will act as the gate
    # A and B are the inputs, they merge into AB to activate the promoter and produce the protein
    A, B, AB, Promoter, Protein = BaseSpecies()
    A + B >> AB [0.5]
    AB + Promoter.inactive >> Promoter.active [0.5]
    Promoter >> Promoter + Protein [lambda promoter: 1 if promoter.active else 0]
    Protein >> Zero [2]

    Promoter(100)
    A(A_conc), B(B_conc)
    MySim = Simulation(A | B | AB | Promoter | Protein)
    MySim.duration = 10
    MySim.level = 0
    MySim.save_data = False
    MySim.plot_data = False
    MySim.run()
    return MySim.results[Protein][0][-1]


# Here is a simple heatmap definition to show to get working
heatmap = []
for a in [0, 25, 50, 75, 100]:
    heatmap_line = []
    for b in [0, 25, 50, 75, 100]:
        output = AND_GATE(a, b)
        heatmap_line.append(output)
    heatmap.append(heatmap_line)

for line in heatmap:
    print(line)

fig = px.imshow(heatmap, x=[0, 25, 50, 75, 100], y=[0, 25, 50, 75, 100], labels=dict(x='b', y='a'))
fig.show()
