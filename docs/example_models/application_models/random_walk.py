import sys, os
import matplotlib.pyplot as plt
from mobspy import *
import time


"""
    This is a geometry-based model
    Here the Meta-Species Mesh is used to represent the positions in a grid
    And their transitions using a CRN
    Here we have a bacteria and a phage that move through a Mesh and once they encounter the bacteria becomes infected
    Although the movement is not realistic, it shows MobsPy geometry capabilities 
    In future work we hope to implement different type movements or grids
"""

start_time = time.time()

Mesh = BaseSpecies()
n = 5
for i in range(n):
    for j in range(n):
        coordinate = 'p_' + str(i) + '_' + str(j)

        if i + 1 < n:
            Mesh.c(coordinate) >> Mesh.c(f'p_{i+1}_{j}')[0.1]
        if i - 1 > -1:
            Mesh.c(coordinate) >> Mesh.c(f'p_{i-1}_{j}')[0.1]
        if j - 1 > -1:
            Mesh.c(coordinate) >> Mesh.c(f'p_{i}_{j-1}')[0.1]
        if j + 1 < n:
            Mesh.c(coordinate) >> Mesh.c(f'p_{i}_{j+1}')[0.1]

Bacteria, Phage = New(Mesh)
Bacteria.not_infected + Phage >> Bacteria.infected [lambda r1, r2: 1000000 if Mesh(r1) == Mesh(r2) else 0]
Bacteria.p_0_0(1)
Phage.c(f'p_{n-1}_{n-1}')(1)

MySim = Simulation(Bacteria | Phage)
MySim.simulation_method = 'stochastic'
MySim.repetitions = 1
MySim.duration = 100
MySim.output_event = True
MySim.plot_data = False
MySim.save_data = False

MySim.run()


def grab_position(species_position, list_x, list_y):
    _, x, y = species_position.split('_')
    list_x.append(int(x))
    list_y.append(int(y))


data = MySim.fres

Bacteria_x = []
Bacteria_y = []
Phage_x = []
Phage_y = []
Meeting_x = []
Meeting_y = []
break_flag = False
for t in range(len(data['Time'])):
    position_flag_bacteria = False
    position_flag_phage = False

    for key in data:

        if key == 'Time' or key == 'Bacteria' or key == 'Phage':
            continue
        species_string = deepcopy(key)

        if data[key][t] == 1:

            split_key = species_string.split('.')

            if split_key[0] == 'Bacteria':
                if split_key[1] == 'infected':
                    grab_position(split_key[-1], Meeting_x, Meeting_y)
                    break_flag = True
                    break

                if not position_flag_bacteria:
                    grab_position(split_key[-1], Bacteria_x, Bacteria_y)
                    position_flag_bacteria = True
                else:
                    raise TypeError('Something went wrong')

            if split_key[0] == 'Phage':
                if not position_flag_phage:
                    grab_position(split_key[-1], Phage_x, Phage_y)
                    position_flag_phage = True
                else:
                    raise TypeError('Something went wrong')

    if break_flag:
        break

plt.plot(Bacteria_x, Bacteria_y, c='b')
plt.plot(Phage_x, Phage_y, c='r')
plt.plot(Meeting_x, Meeting_y, c='y')
plt.show()
