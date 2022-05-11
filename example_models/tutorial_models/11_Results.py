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

# For the data plot, MobsPy has a default plot
# The default plot changes for both stochastic and deterministic
MySim.plot_data = False

# Run the sim
MySim.run()

# For the deterministic plot call
# If no species are given, MobsPy will plot them all
MySim.plot_stochastic(A, C)

# If save_data is true MobsPy will generate a JSON file containing the data
# Otherwise it can be accessed directly from the simulation object
# The simulation results contain the data, parameters and mapping used in the model
# The data is separated by species name as keys
# And stored with all the runs in the runs key - to access them use natural numbers
print(MySim.results['data']['A']['runs'][0])

"""
    Plotting with queries
"""

# If the meta-species has several different species inside - MobsPy will plot the sum as default
# If only a certain characteristic wishes to be plot just query for it inside the plot
# It will sum over all species with that characteristic

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

