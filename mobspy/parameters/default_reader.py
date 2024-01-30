"""
    This model is responsible for storing MobsPy default parameters
"""
from mobspy.parameter_scripts.parameter_reader import read_json


# This was just created to avoid potential directory compatibilities
def get_default_parameters():
    default_parameters = {
        "__comment_1": "Model parameters - Repetitions only for stochastic",

        "volume": 1,
        "repetitions": 1,
        "level": 3,
        "rate_type": None,

        "__comment_2": "basiCO parameters",

        "simulation_method": "deterministic",
        "method": None,
        "start_time": 0,
        "duration": 60,
        "r_tol": 1e-8,
        "a_tol": 1e-10,

        "__comment_3": "Core number",

        "jobs": -1,

        "__comment_4": "OUTPUT",

        "output_dir": "outputs/",
        "output_file": None,
        "output_event": False,
        "unit_x": None,
        "unit_y": None,
        "skip_expression_check": False,
        "output_concentration": True,
        "save_data": False,
        "plot_data": True,
        "plot_type": None,

        "_continuous_simulation": False,
        "_end_condition": None,
        "_with_event": False
    }
    return default_parameters
