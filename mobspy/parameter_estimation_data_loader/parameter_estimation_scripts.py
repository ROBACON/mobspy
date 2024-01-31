import mobspy.simulation_logging.log_scripts as simlog
import basico
from pandas import DataFrame


def basiCO_parameter_estimation(simulation_object, parameters_to_estimate,
                                experimental_data=None, bound=(1e-3, 1e3),
                                method='Evolution Strategy (SRES)'):
    """
        This function fits a MobsPy parameter in a model to experimental data

        :param simulation_object: (Simulation Object) - MobsPy Simulation Object
        :param parameters_to_estimate: (list of MobsPy parameters, or list of str) List of MobsPy parameters
        to estimate
        :param experimental_data: (list of pandas Dataframe, or pandas Dataframe) - experimental data to fit the
        parameter with
        :param bound: (list, tuple) - list of two elements with a lower and upper bound of the parameters values
        :param method: (str) - method to be used by basiCO optimisation - The options are Random Search,
        Simulated Annealing, Differential Evolution, Scatter Search, Genetic Algorithm, Evolutionary Programming,
        Genetic Algorithm SR, Evolution Strategy (SRES), Particle Swarm
    """

    # Check of valid inputs
    if simulation_object.experimental_data is not None:
        experimental_data = simulation_object.experimental_data.return_pandas()[0]
    if experimental_data is None:
        simlog.error('No experimental data found in the simulation object or as argument of the '
                     'basiCO_parameter_estimation function')

    if type(experimental_data) != list and type(experimental_data) != set and type(experimental_data) != tuple \
            and not isinstance(experimental_data, DataFrame):
        simlog.error('Experimental for basiCO estimation must be a list of pandas dataframes or a pandas dataframe')

    if type(parameters_to_estimate) != list and type(parameters_to_estimate) != set \
            and type(parameters_to_estimate) != tuple:
        simlog.error('The parameter that will be estimated must be inside a list, set or tuple')

    # Conversion of MobsPy parameters to str
    converted_parameters = []
    for par in parameters_to_estimate:
        converted_parameters.append(str(par))
    parameters_to_estimate = converted_parameters

    sbml_list = simulation_object.generate_sbml()
    if len(sbml_list) > 1:
        simlog.error("BasiCO optimization does not support composite simulation optimization")

    sbml_str = sbml_list[0]
    basico.model_io.load_model_from_string(sbml_str)

    fit_list = []
    for par in parameters_to_estimate:
        basico_reaction_dict = find_parameters_in_basico_dataframe(basico.get_reaction_parameters(), par).to_dict()
        basico_parameter_name = list(basico_reaction_dict['reaction'].keys())[0]

        fit_dictionary = {'name': basico_parameter_name, 'lower': bound[0], 'upper': bound[1]}
        fit_list.append(fit_dictionary)

    if type(experimental_data) == list or type(experimental_data) == set or type(experimental_data) == tuple:
        for i, exp in enumerate(experimental_data):
            basico.add_experiment('exp' + str(i), exp)
    else:
        basico.add_experiment('exp1', experimental_data)

    basico.set_fit_parameters(fit_list)
    basico_results = basico.run_parameter_estimation(method=method)

    # WHY IS EVERYTHING PANDAS IN BASICO??????????? - I don't get paid enough for this
    results = {}
    for key, sol in basico_results.to_dict()['sol'].items():
        parameter_name = str(key).replace('Values[', '')
        parameter_name = str(parameter_name).replace(']', '')
        results[parameter_name] = sol
    print(results)

    exit()


def find_parameters_in_basico_dataframe(basico_reactions_df, mobspy_parameter_name):
    """
        Finds the corresponding mobspy parameter in the basico reactions parameters dataframe

        :param basico_reactions_df: (Dataframe) basiCO reactions df
        :param mobspy_parameter_name: (str) name of the mobspy parameter
    """

    def has_mobspy_parameter_in_name(basico_reaction_name, mobspy_parameter_name):
        # Change this function based on your specific logic for finding common substrings
        parameter_name = basico_reaction_name.split('.')[1]
        return True if parameter_name == mobspy_parameter_name else False

    def find_common_substrings(df):
        return has_mobspy_parameter_in_name(df, mobspy_parameter_name)

    return basico_reactions_df[basico_reactions_df.index.to_frame()['name'].apply(find_common_substrings)]


def python_parameter_estimation():
    # Work in progress
    pass
