from mobspy import *

'''
    Here we will how MobsPy handles results
    Both as a plotting aspect and result accessing aspect
'''

# Firstly define the model
A, B, C = BaseSpecies(3)
A(200) + B(100) >> C [0.001]
MySim = Simulation(A | B | C)
MySim.save_data = False

# Set to stochastic
MySim.simulation_method = 'stochastic'

# MobsPy usually has two default plots
# One for stochastic and one for deterministic
# The default plots will be called accordingly to the simulation_method
# after the simulation
# Unless plot_data is set to false
MySim.plot_data = False

# Run the sim
MySim.run()

# If the default plot was not used one can always call the plot_stochastic or plot_deterministic
# from the resulting simulation object
# These function take either species as arguments or nothing at all
# If no argument is passed all species will be plotted
# If some species are passed as argument, they will be the only ones to be plotted
# For the stochastic_plot the default plot consists of two figures for each species
# One for the runs and the other for the average value with standard deviation
MySim.plot_stochastic(A, C)

# If save_data is true MobsPy will generate a JSON file containing the data
# Otherwise it can be accessed from the simulation object
# The simulation results contain the data, parameters and mappings used in the model
# The data is separated by species name
# And stored with all the runs under the runs key
# For example to recover the first run for the A species:
print(MySim.results['data']['A']['runs'][0])

"""
    Plotting with queries
"""


# Plotting functions can also query
# For instance in the model bellow we are interested in the evolution of a species with a certain characteristic
# We want to observe all the values of C with the characteristic b1
# So we just pass C.b1 as an argument to the plotting function

A, B = BaseSpecies(2)
A.a1, A.a2
B.b1, B.b2
C = A*B
C.b1.a1(100), C.b1.a2(100)
C.b1 >> C.b2 [0.1]
MySim = Simulation(C)
MySim.save_data = False
MySim.plot_data = False
MySim.run()

MySim.plot_deterministic(C.b1)

