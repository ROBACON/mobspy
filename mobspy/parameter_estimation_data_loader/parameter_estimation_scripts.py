import mobspy.simulation_logging.log_scripts as simlog
from mobspy.import_manager.lazy_import_class import LazyImporter as ipm_LazyImporter
basico = ipm_LazyImporter('basico')
mobspy_basico_patch = ipm_LazyImporter('mobspy.patch_scripts.basico_task_parametrization')
from pandas import DataFrame


def basiCO_parameter_estimation(simulation_object, parameters_to_estimate,
                                experimental_data=None, bound=None,
                                method='Evolution Strategy (SRES)', verbose=True,
                                change_parameter_values=True):
    """
        This function fits a MobsPy parameter in a model to experimental data.
        Bound specifies the search range.
        If the bound is None it takes as bounds /1000 and *1000 the original value of the parameter
        The function updates the value of the parameter when it is done

        :param simulation_object: (Simulation Object) - MobsPy Simulation Object
        :param parameters_to_estimate: (list of MobsPy parameters, or list of str) List of MobsPy parameters
        to estimate
        :param experimental_data: (list of pandas Dataframe, or pandas Dataframe) - experimental data to fit the
        parameter with
        :param bound: (list, tuple, dict) - list of two elements with a lower and upper bound of the parameters values,
        or dictionary with the parameter name and a two element for that specific parameter
        :param method: (str) - method to be used by basiCO optimisation - The options are Random Search,
        Simulated Annealing, Differential Evolution, Scatter Search, Genetic Algorithm, Evolutionary Programming,
        Genetic Algorithm SR, Evolution Strategy (SRES), Particle Swarm
        :param verbose: (bool) print the results after finishing or not
        :param change_parameter_values: (bool) change/convert parameter values when possible
    """

    # Check inputs in order ###############################################################################
    # Parameters ok?
    if type(parameters_to_estimate) != list and type(parameters_to_estimate) != set \
            and type(parameters_to_estimate) != tuple:
        simlog.error('The parameter that will be estimated must be inside a list, set or tuple')

    converted_parameters = []
    # If bound is None, we set the auto-bound when checking the parameters
    if bound is None:
        bound = {}
        flag_auto_set = True
    else:
        flag_auto_set = False

    # Conversion of MobsPy parameters to str, construction of bound
    original_parameters = parameters_to_estimate
    for par in parameters_to_estimate:
        if flag_auto_set:
            bound[str(par)] = [par.value/1000, par.value*1000]

        converted_parameters.append(str(par))
    parameters_to_estimate = converted_parameters

    # Simulation Object ok?
    if simulation_object.experimental_data is not None:
        experimental_data = simulation_object.experimental_data.return_pandas()[0]
    if experimental_data is None:
        simlog.error('No experimental data found in the simulation object or as argument of the '
                     'basiCO_parameter_estimation function')

    # Experimental data ok?
    if type(experimental_data) != list and type(experimental_data) != set and type(experimental_data) != tuple \
            and not isinstance(experimental_data, DataFrame):
        simlog.error('Experimental for basiCO estimation must be a list of pandas dataframes or a pandas dataframe')

    # Bound ok?
    new_bound = {}
    if type(bound) == dict:
        for key in bound:
            new_bound[str(key)] = bound[key]
        bound = new_bound

        for par in parameters_to_estimate:
            if par not in bound:
                simlog.error('If a dictionary is used for the bounds, all parameters range for estimation '
                             'must be specified. Make sure the dictionary keys contains all parameters and a list with '
                             'upper and lower bound value for the each parameter is given as the items ')
    else:
        try:
            if len(bound) != 2:
                simlog.error('The bound argument must be a list with the lower and upper bound of all parameters')
        except:
            simlog.error('The bound argument must be a list with the lower and upper bound of all parameters')
    ########################################################################################################

    sbml_list = simulation_object.generate_sbml()
    if len(sbml_list) > 1:
        simlog.error("BasiCO optimization does not support composite simulation optimization")

    sbml_str = sbml_list[0]
    model = basico.model_io.load_model_from_string(sbml_str)

    fit_list = []
    for par in parameters_to_estimate:

        basico_reaction_dict = find_parameters_in_basico_dataframe(basico.get_reaction_parameters(), par).to_dict()
        try:
            basico_parameter_name = list(basico_reaction_dict['reaction'].keys())[0]
        except IndexError:
            simlog.error(f'Parameter {par} was not found in the Simulation model. \n '
                         f'Please make sure that any of the meta-species used to construct the simulator use '
                         f' the parameter in one of their reactions.')

        if bound is None:
            pass

        if type(bound) == dict:
            fit_dictionary = {'name': basico_parameter_name, 'lower': bound[par][0], 'upper': bound[par][1]}
        else:
            fit_dictionary = {'name': basico_parameter_name, 'lower': bound[0], 'upper': bound[1]}
        fit_list.append(fit_dictionary)

    if type(experimental_data) == list or type(experimental_data) == set or type(experimental_data) == tuple:
        for i, exp in enumerate(experimental_data):
            mobspy_basico_patch.add_experiment('exp' + str(i), exp, model=model)
    else:
        mobspy_basico_patch.add_experiment('exp1', experimental_data, model=model)

    basico.set_fit_parameters(fit_list, model=model)
    basico_results = basico.run_parameter_estimation(model=model, method=method)

    # WHY IS EVERYTHING PANDAS IN BASICO??????????? - I don't get paid enough for this
    results = {}
    for key, sol in basico_results.to_dict()['sol'].items():
        parameter_name = str(key).replace('Values[', '')
        parameter_name = str(parameter_name).replace(']', '')
        results[parameter_name] = sol

    if change_parameter_values:
        for p in original_parameters:
            if p.has_units():
                p.set_value(results[str(p)]).convert_to_original_unit()
                results[str(p)] = p.value
            else:
                p.set_value(results[str(p)])

    if verbose:
        print('Parameter estimation complete. The results follow: ')
        for key in results:
            print(key, results[key])

    return results


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
