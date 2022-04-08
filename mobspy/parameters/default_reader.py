from mobspy.parameter_scripts.parameter_reader import read_json


# This was just created to avoid potential directory compatibilities
def get_default_parameters():
    default_parameters = {
        "__comment_1": "Model parameters - Repetitions only for stochastic",

        "volume": 1,
        "repetitions": 3,
        "level": 3,

        "__comment_2": "basiCO parameters",

        "simulation_method": "deterministic",
        "start_time": 0,
        "duration": 60,
        "r_tol": 1e-8,
        "a_tol": 1e-10,

        "__comment_3": "OUTPUT",

        "output_dir": "outputs/",
        "output_file": None,
        "output_event": False,
        "unit_x":None,
        "unit_y":None,
        "output_concentration":False,
        "save_data": True,
        "plot_data": True

    }
    return default_parameters
