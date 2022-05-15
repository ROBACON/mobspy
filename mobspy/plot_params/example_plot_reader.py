"""
    This model is responsible for storing example plotting parameters
    It's used by the simulation to check which parameters exist or not
"""

def get_example_plot_parameters():
    example_parameters = {
        "output_dir": "",
        "logscale": ["X", "Y"],
        "xlim": [0, 1],
        "ylim": [0, 1e3],
        "figsize": [1.5, 1.5],
        "no_title": True,
        "time_series": [0],
        'linestyle':'-',
        'linewidth':0.1,
        'label':'test_label',
        'runs':[0],
        'fill_between':False,
        'frameon':False,
        'title':'Test_Plot',
        'annotations':['Test_note', 0, 0],
        'pad':1,
        'dpi':3,
        "plots":False,
        "species_to_plot": None,

        "figures": [{
            "plots": [{
                "species_to_plot": ["A"]
            }]
        }],

        # These parameters are not supported by plot_raw
        # Their only function is to alter the plot labels so please use that parameter
        # They are only listed here to avoid printing out errors during parameter checks
        "unit_x": '',
        "unit_y": '',
        "output_concentration": False,
        'simulation_method':'stochastic'
    }
    return example_parameters

