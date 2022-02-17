from parameter_scripts.parameter_reader import read_json


# This was just created to avoid potential directory compatibilities
def get_default_parameters():
    default_parameters = {
        "__comment_1": "Model parameters - Repetitions only for stochastic",

        "volume_ml": 1,
        "repetitions": 3,

        "__comment_2": "basiCO parameters",

        "simulation_method": "deterministic",
        "start_time": 0,
        "duration": None,
        "r_tol": 1e-8,
        "a_tol": 1e-10,

        "__comment_3": "OUTPUT",

        "output_dir": "outputs/",
        "output_file": None,
        "save_data": True,
        "plot": True

    }
    return default_parameters
