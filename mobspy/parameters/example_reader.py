"""
    This model is responsible for storing example parameters
    It's used by the simulation to check which parameters exist or not.
    If a new parameter is added, they key must be added here as an example.
    Otherwise the simulation object will throw a compilation error, telling the user that the parameter is not
    supported,
"""


def get_example_parameters():
    example_parameters = {
        "__comment_1": "Model parameters - Repetitions only for stochastic",

        "volume": 1,
        "repetitions": 3,
        "level": 0,

        "__comment_2": "basiCO parameters",

        "simulation_method": "stochastic",
        "method": None,
        "start_time": 0,
        "duration": None,
        "r_tol": 1e-8,
        "a_tol": 1e-10,
        "seeds": [1, 2, 3],
        "step_size": 1,

        "__comment_3": "OUTPUT",

        "unit_x": "year",
        "unit_y": "nanomolar",
        "output_concentration": True,
        "output_dir": "",
        "output_event": False,
        "output_file": None,
        "save_data": False,
        "plot_data": True

    }
    return example_parameters
