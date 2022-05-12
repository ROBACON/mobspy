from mobspy import *

A, B, C = BaseSpecies(3)
A(200) + B(100) >> C[0.005]
MySim = Simulation(A | B | C)
MySim.save_data = False
MySim.plot_data = False

# The hierarchical plot structure was designed to facilitate the plotting of multiple
# figures and curves within a single plot window
#
# In this plotting structure there is a figures key that accepts a list.
# For each element in this list, the hierarchical plotting structure will create a new figure
# Each figure contains a plots list with the number of curves to be plotted inside the figure
#
# Now the hierarchical plot structure works in the following manner:
#   It follows the order:
#   => plot-wise, then figure-wise and finally global-wise
#   The configuration parameters can be placed inside the plots, figures or as a key to the dictionary itself
#   When plotting MobsPy will check to see if the parameter was defined for the curve
#   If it was defined for the curve the parameter is used for plotting
#   If it is not defined it checks the figure
#   Then it repeats the same process for figure-wise before moving to global-wise
#
# Also the plotting structure allows one to assign some parameters by species
# To do so just create a dictionary where the key is the species name as can be seen bellow
# Based on this, if you wish to use MobsPy plotting structure please avoid naming your species
# with plotting parameters.
# A viable strategy is to always capitalize species names as no parameter is capitalized
#
# Both plot_stochastic and plot_deterministic can be configured using this structure
# Plot_raw can only be configured through this structure and does not accept the dot-based configuration
example_parameters = {

    'A': {'color': 'red',  'label': 'A'},
    'B': {'color': 'blue',  'label': 'B'},
    'C': {'color': 'green', 'label': 'C'},
    "species_to_plot":["A"],

    'figures':[{
            "plots":[
                {"species_to_plot": ["B"]},
                {}
        ]},{
            "plots":[
                {"species_to_plot": ["C"]},
                {}
        ]}
    ]
}

MySim.run()
MySim.plot_raw(example_parameters)
