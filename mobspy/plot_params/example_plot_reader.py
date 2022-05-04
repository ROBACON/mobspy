
def get_example_plot_parameters():
    example_parameters = {
        "output_dir": "",
        "logscale": ["X", "Y"],
        "xlim": [0, 1],
        "ylim": [0, 1e3],
        "figsize": [1.5, 1.5],
        "no_title": True,

        "figures": [{
            "plots": [{
                "species_to_plot": ["A"]
            }]
        }]
    }
    return example_parameters