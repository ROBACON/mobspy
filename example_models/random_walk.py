import sys, os
import matplotlib.pyplot as plt
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *

Mesh = BaseSpecies(1)
n = 10
for i in range(n):
    for j in range(n):
        coordinate = 'p_' + str(i) + '_' + str(j)

        if i + 1 < n:
            Mesh.c(coordinate) >> Mesh.c(f'p_{i+1}_{j}') [0.1]
        if i - 1 > -1:
            Mesh.c(coordinate) >> Mesh.c(f'p_{i-1}_{j}') [0.1]
        if j - 1 > -1:
            Mesh.c(coordinate) >> Mesh.c(f'p_{i}_{j-1}') [0.1]
        if j + 1 < n:
            Mesh.c(coordinate) >> Mesh.c(f'p_{i}_{j+1}') [0.1]

Bacteria = Mesh*New
Phage = Mesh*New
Bacteria.not_infected + Phage >> Bacteria.infected [lambda r1, r2: 1000000 if Mesh(r1) == Mesh(r2) else 0]
Bacteria.p_0_0(1)
Phage.p_9_9(1)

MySim = Simulation(Bacteria | Phage)
MySim.simulation_method = 'stochastic'
MySim.repetitions = 1
MySim.duration = 100
MySim.output_event = True
MySim.plot_data = False
MySim.save_data = False

MySim.run()


def grab_position(species_position, list_x, list_y):
    _, x, y = Mesh(species_position).split('_')
    list_x.append(int(x))
    list_y.append(int(y))


data = MySim.results['data']
print(data['Time'])

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

        if data[key]['runs'][0][t] == 1:

            if species_string.split('.')[0] == 'Bacteria':
                if Bacteria(key) == 'infected':
                    grab_position(key, Meeting_x, Meeting_y)
                    break_flag = True
                    break

                if not position_flag_bacteria:
                    grab_position(key, Bacteria_x, Bacteria_y)
                    position_flag_bacteria = True
                else:
                    raise TypeError('Something went wrong')

            if species_string.split('.')[0] == 'Phage':
                if not position_flag_phage:
                    grab_position(key, Phage_x, Phage_y)
                    position_flag_phage = True
                else:
                    raise TypeError('Something went wrong')

    if break_flag:
        break

plt.scatter(Bacteria_x, Bacteria_y, marker='o', c='b')
plt.scatter(Phage_x, Phage_y, marker='o', c='r')
plt.scatter(Meeting_x, Meeting_y, marker='o', c='y')
plt.show()
