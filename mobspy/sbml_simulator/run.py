# import basico
# Lazy import basico
from mobspy.import_manager.lazy_import_class import LazyImporter as ipm_LazyImporter
basico = ipm_LazyImporter('basico')
from joblib import Parallel, delayed
import mobspy.simulation_logging.log_scripts as simlog
import mobspy.sbml_simulator.builder as sbml_builder
from copy import deepcopy


# def simulate(list_of_params, models)
def simulate(jobs, list_of_params, models):
    """
        This function coordinates the simulation by calling the necessary jobs
        In the future we hope to implement parallel cluster computing compatibility

        :param list_of_params: (dict) simulation parameters from the text file
        :param models: (list) [{'species_for_sbml':, 'parameters_for_sbml':, 'reactions_for_sbml':,
                            'events_for_sbml':, 'species_not_mapped':, 'mappings':}]
        :return: data (dict) = dictionary containing the resulting data from simulation
    """

    # Run in parallel or sequentially
    # If nothing is specified just run it in parallel
    data = job_execution(list_of_params, models, jobs)

    return data


def job_execution(params, models, jobs):
    """
        This is defined for parallelism purposes
        Uses multiple cores from the processor to execute stochastic simulations

        :param params: (dict) = simulation parameters
        :param models: (list) [{'species_for_sbml':, 'parameters_for_sbml':, 'reactions_for_sbml':,
                            'events_for_sbml':, 'species_not_mapped':, 'mappings':}]
        :param jobs: (int) = number of cores to use, -1 for all available

        :return: parallel_data data of all the individual simulations executed in parallel
    """

    def __single_run(packed):
        i = packed

        added_data = {}
        # This is only need to avoid git errors
        reformatted_data = {}
        for j, (sim_par, model) in enumerate(zip(params, models)):

            # Generate SBML here
            if j > 0:
                sbml_str = __sbml_new_initial_values(reformatted_data, model, sim_par, new_model=True)
            else:
                sbml_str = __sbml_new_initial_values({}, model, sim_par)

            end_condition_not_satisfied = True
            if sim_par['_continuous_simulation']:
                duration = float(sim_par['initial_conditional_duration'])
            else:
                duration = float(sim_par['duration'])

            while end_condition_not_satisfied:
                basico_model = basico.model_io.load_model_from_string(sbml_str)
                data = __run_time_course(basico_model, duration, sim_par, i)

                reformatted_data = reformat_time_series(data)

                if sim_par['_continuous_simulation']:
                    if reformatted_data['End_Flag_MetaSpecies'][-1] > 0:
                        reformatted_data = __filter_condition_event_time_data(reformatted_data)
                        end_condition_not_satisfied = False
                    else:
                        sbml_str = __sbml_new_initial_values(reformatted_data, model, sim_par)
                        duration = 2 * duration
                else:
                    end_condition_not_satisfied = False

                reformatted_data = __remap_species(reformatted_data, model['mappings'], model['species_for_sbml'])
                added_data = __add_simulations_data(added_data, reformatted_data)

        return added_data

    parallel_data = Parallel(n_jobs=jobs)(delayed(__single_run)(i) for i in range(params[0]['repetitions']))

    if not parallel_data:
        simlog.error("Error: The parallel model has not produced an output." +
                     "Try addding ('sequential': True) to parameters")

    return parallel_data


def __run_time_course(basico_model, duration, params, index):
    """
        This function prepares the parameters for the basiCO.run_time_course

        :param basico_model: basico model to be run in the simulation
        :param params: (dict) simulation parameters
        :param index: (int) current index of the run
        :param duration: (int, float) simulation duration in seconds
    """
    params["simulation_method"] = params["simulation_method"].lower()
    if (params['_with_event'] or params['_continuous_simulation']) and params["simulation_method"] == 'stochastic':
        params["simulation_method"] = 'directmethod'

    kargs = {'model': basico_model,
             'method': params["simulation_method"],
             'start_time': params["start_time"],
             'r_tol': params["r_tol"],
             'a_tol': params["a_tol"],
             'output_event': params['output_event']}

    if 'seeds' in params:
        kargs['use_seed'] = True
        kargs['seed'] = params['seeds'][index]

    if 'step_size' in params:
        kargs['automatic'] = False
        kargs['step_number'] = int(params['duration'] / params['step_size'])

    return basico.run_time_course(duration, **kargs)


def reformat_time_series(data):
    """
        Transforms the _dot_ in from the results into .

        :param data: (pd.dataframe) simulation data results from BasiCO
    """
    data_dict = {'Time': data.index.tolist()}

    for key in data:
        data_dict[key.replace('_dot_', '.')] = list(data[key])

    return data_dict


def __filter_condition_event_time_data(data):
    """
        This function filters all the data until the event was triggered. It is used for the conditional duration
        simulations

        :param data: (dict) single run data with species as keys and values
    """
    new_data = {}

    for i, e in enumerate(data['End_Flag_MetaSpecies']):
        if e == 1:
            stop_index = i
            break

    for key in data:
        new_data[key] = data[key][:stop_index + 1]

    return new_data


def __sbml_new_initial_values(data, model, sim_para, new_model=False):
    """
        This function adjusts the new sbml counts to the end of the previous concatenation simulation or conditional
        simulation

        :param data: (dict) single run data with species as keys and values
        :param model: (dict) model to be used in the next simulation
        :param model: (sim_par) parameters to be used in the next simulation
        :param new_model: (bool) is a new model being used or the model is the same as the old one
    """
    species_for_sbml = model['species_for_sbml']

    check_list = ["stochastic", "directmethod"]
    for key in data:
        sbml_key = key.replace('.', '_dot_')
        if sbml_key not in species_for_sbml.keys():
            continue

        if key == 'Time':
            continue
        try:
            # Case of species set
            if sim_para["simulation_method"].lower() in check_list:
                species_for_sbml[sbml_key] = int(list(data[key])[-1])
            else:
                species_for_sbml[sbml_key] = list(data[key])[-1]
        except KeyError:
            pass

    if new_model:
        try:
            species_for_sbml['End_Flag_MetaSpecies'] = 0
        except KeyError:
            pass

    return sbml_builder.build(species_for_sbml, model['parameters_for_sbml'],
                              model['reactions_for_sbml'], model['events_for_sbml'],
                              model['assignments_for_sbml'])


def __add_simulations_data(added_data, reformatted_data):
    """
        Adds the data from the new executed simulation to the data stored so far

        :param added_data: (dict) added data so far in the simulation (key: species string - value: run)
        :param reformatted_data: (dict) data return from the simulation after being formatted using the reformat
        data function
    """
    if added_data != {}:
        time_to_add = added_data['Time'][-1]
    else:
        time_to_add = 0
    new_data = {}
    already_added_keys = set()

    for key in added_data:
        new_data[key] = added_data[key]

    for i, time in enumerate(reformatted_data['Time']):
        # Remove the repeated initial value from following simulations.
        # The final value of summed simulations is repeated
        if time == 0 and added_data != {}:
            for key in reformatted_data:
                reformatted_data[key].pop(0)

        # This is for simulations with only one value. The system did not change at all
        try:
            reformatted_data['Time'][i] = reformatted_data['Time'][i] + time_to_add
        except IndexError:
            pass

    for key in added_data:
        if key == 'Time':
            continue
        try:
            new_data[key] = added_data[key] + reformatted_data[key]
            already_added_keys.add(key)
        except KeyError:
            dummy = added_data[key][-1]
            new_data[key] += [dummy for _ in reformatted_data['Time']]

    for key in reformatted_data:
        if key == 'Time':
            continue

        if key in already_added_keys:
            continue
        else:
            if time_to_add != 0:
                new_data[key] = [0 for _ in added_data['Time']]
                new_data[key] = new_data[key] + reformatted_data[key]
            else:
                new_data[key] = reformatted_data[key]

    if time_to_add != 0:
        new_data['Time'] = added_data['Time'] + reformatted_data['Time']
    else:
        new_data['Time'] = reformatted_data['Time']

    return new_data


def __remap_species(data, mapping, species_not_mapped):
    """
        Takes the simulated species (data) and add ones defined by
        mapping (mapping).
        By default, a mapping is a list of species in which case the
        sum is taken: ['a', 'b', ...]
        It also adds species that have been removed from the basico simulation

        Parameters:
        :param data: (dict) Data in MobsPy format
        :param mapping: (dict) = Mappings between species and meta-species
        :param species_not_mapped: (list of str) Species not mapped by BasiCO
    """

    mapped_data = {'Time': data['Time']}
    T = range(len(data['Time']))

    # copy over all unmapped ones
    for k in data.keys():
        mapped_data[k] = data[k]

    dot_species_not_mapped = {}
    for key in species_not_mapped:
        dot_species_not_mapped[key.replace('_dot_', '.')] = species_not_mapped[key]

    # 1st pass with sum mappings
    for group in mapping.keys():

        the_mapping = mapping[group]
        mapped_data[group] = {'runs': []}

        try:
            # check if is a list -> sum
            if type(the_mapping) is list:

                this_run = []
                runs_not_returned_by_basico = {}
                for t in T:
                    mapping_sum = 0
                    for spe in the_mapping:
                        try:
                            mapping_sum = mapping_sum + data[spe][t]
                        except KeyError:
                            mapping_sum = mapping_sum + dot_species_not_mapped[spe]
                            try:
                                runs_not_returned_by_basico[spe] += [dot_species_not_mapped[spe]]
                            except KeyError:
                                runs_not_returned_by_basico[spe] = [dot_species_not_mapped[spe]]
                    this_run.append(mapping_sum)
                for spe in runs_not_returned_by_basico:
                    try:
                        mapped_data[spe] = runs_not_returned_by_basico[spe]
                    except KeyError:
                        mapped_data[spe] = [runs_not_returned_by_basico[spe]]
                mapped_data[group] = this_run

        except IndexError:
            simlog.error(f'run: remap_species: error when remapping "{the_mapping}".' +
                         'Possible fix: All runs must have the same time')

        except TypeError:
            simlog.warning(f'Copasi removes A >> A species from reaction calculations and does not provide an output')
            simlog.warning(f'Please check the output data to see if this is the problem')
            for key in data:
                print(key, data[key])
            exit(1)

    return mapped_data
