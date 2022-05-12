from mobspy import *

if __name__ == '__main__':

    A, B, C = BaseSpecies(3)
    A(200) + B(100) >> C[0.005]
    MySim = Simulation(A | B | C)
    MySim.save_data = False
    MySim.plot_data = False
    MySim.unit_x = 'minutes'

    MySim.run()
    MySim.plot_deterministic(A)
    exit()

    example_parameters = {

        'A': {'color': 'red', 'label': 'A'},
        'B': {'color': 'blue', 'label': 'B'},
        'C': {'color': 'green', 'label': 'C'},
        "species_to_plot": ["A"],

        'figures': [{
            "plots": [
                {"species_to_plot": ["B"]},
                {}
            ]}, {
            "plots": [
                {"species_to_plot": ["C"]},
                {}
            ]}
        ]
    }

    # MySim.run()
    # MySim.plot_raw(example_parameters)