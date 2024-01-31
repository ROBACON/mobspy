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
        'runs': [0],
        'fill_between': False,
        'frameon': False,

        'xlabel': 'Something',
        'ylabel': 'Something Else',
        'xlabel_fontsize': 14,
        'ylabel_fontsize': 14,
        'title': 'Test_Plot',
        'title_fontsize': 16,
        'suptitle': 'Test_Sup_Title',
        'suptitle_fontsize': 18,

        'annotations': ['Test_note', 0, 0],
        'pad': 1,
        'dpi': 3,
        "plots": False,
        "species_to_plot": None,
        "save_to": 'example.png',

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
        'simulation_method': 'stochastic',
        "time_filter": [0, 100000]
    }
    return example_parameters

