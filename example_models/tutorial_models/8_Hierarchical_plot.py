from mobspy import *

A, B, C = BaseSpecies(3)
A(200) + B(100) >> C[0.005]
MySim = Simulation(A | B | C)
MySim.save_data = False
MySim.plot_data = False

# The hierarchical plot structure looks for parameters in the following order
# First plot-wise, then figure-wise and finally global-wise
# The override rule is plot-wise > figure-wise > global-wise
# If the parameter is not found at a local value it searches above
# Finally, plot raw only accepts dictionaries and JSON
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
