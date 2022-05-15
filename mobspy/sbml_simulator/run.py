import basico
from joblib import Parallel, delayed
import mobspy.simulation_logging.log_scripts as simlog


def simulate(sbml_str, params, mappings, species_not_mapped):
    """
        This function coordinates the simulation by calling the necessary jobs
        In the future we hope to implement parallel cluster computing compatibility

        Parameters:
            sbml_str (str) = SBML str from a model for running the Copasi simulator
            mappings (dict) = Mappings that sum the species into the meta-species
            params (dict) = simulation parameters from the text file
            species_for_sbml (dict) = This is needed because basiCO does not return species such that dS/dt = 0

        Returns:
            (dict) Data for MobsPy simulation object
    """

    # Run in parallel or sequentially
    # If nothing is specified just run it in parallel
    try:
        if params["jobs"] == 1:
            simlog.debug("Running simulation sequentially")
            jobs = params["jobs"]
        else:
            simlog.debug("Running simulation in parallel")
            jobs = params["jobs"]
    except KeyError:
        simlog.debug("Running simulation in parallel")
        jobs = -1

    # TODO: If cluster compatibility is added I suggest here
    data = job_execution(sbml_str, params, jobs)
    data = remap_species(data, mappings, params, species_not_mapped)

    simlog.debug("Simulation is Over")
    return {'data': data, 'params': params, 'mappings': mappings}


def job_execution(sbml_str, params, jobs):
    """
        This is defined for parallelism purposes
        Uses multiple cores from the processor to execute stochastic simulations
        
        Parameters:
            params (dict) = simulation parameters
            jobs (int) = number of cores to use, -1 for all available
            sbml_str (str) = model SBML string
    """
    def __single_run(packed):
        sbml_str, i = packed

        basico.model_io.load_model_from_string(sbml_str)
        data = __run_time_course(params['duration'], params, i)

        reformated_data = reformat_time_series(data)

        return reformated_data

    parallel_data = Parallel(n_jobs=jobs)(delayed(__single_run)((sbml_str, i)) for i in range(params['repetitions']))

    if not parallel_data:
        simlog.error("Error: The parallel model has not produced an output." +
                     "Try addding ('sequential': True) to parameters")

    # We always call merge to keep the data in the format we want
    merged_data = merge(params, parallel_data)

    return merged_data


def __run_time_course(duration, params, index):
    """
        This function prepares the parameters for the basiCO.run_time_course

        Parameters:
            params (dict) = simulation parameters
            index (int) = current index of the run
            duration (int, float) = simulation duration in seconds
    """
    kargs = {'method': params["simulation_method"].lower(),
             'start_time': params["start_time"],
             'r_tol': params["r_tol"],
             'a_tol': params["a_tol"],
             'output_event': params['output_event']}

    if 'seeds' in params:
        kargs['use_seed'] = True
        kargs['seed'] = params['seeds'][index]

    if 'step_size' in params:
        kargs['automatic'] = False
        kargs['step_number'] = int(params['duration']/params['step_size'])

    return basico.run_time_course(duration, **kargs)


def reformat_time_series(data):
    """
        Transforms the _dot_ in from the results into .

        Parameters:
            data (pd.dataframe) = simulation data results from BasiCO
    """
    data_dict = {'Time': data.index.tolist()}

    for key in data:
        data_dict[key.replace('_dot_', '.')] = list(data[key])

    return data_dict


def merge(params, data):
    """
    a basic merge of the data with possible resampling in time to save less data

        Parameters
            params (dict) =  dict of parameters
            data (pd.dataframe) =  Array of single runs: data[0], data[1], ...

        Returns:
            merged_data (dict) = Data merged in a dictionary with all time series from all experiments
            separated by species and mapping
    """

    for i in range(1, len(data)):
        assert (data[i - 1]['Time'] == data[i]['Time'])  # should be set at exact same times

    merged_data = {'Time': data[0]['Time']}
    for key in data[0].keys():
        if key not in ['Time']:
            merged_data[key] = {'runs': [data[i][key]
                                         for i in range(params["repetitions"])]}

    return merged_data


def remap_species(data, mapping, params, species_not_mapped):
    """
        Takes the simulated species (data) and add ones defined by
        mapping (mapping).

        By default, a mapping is a list of species in which case the
        sum is taken: ['a', 'b', ...]

        Parameters:
            data (dict) = Data in MobsPy format
            mapping (dict) = Mappings between species and meta-species
            params (dict) = Simulation parameters
            species_not_mapped (list of str) = Species not mapped by BasiCO (removed due to the absence of reactions)
    """

    mapped_data = {'Time': data['Time']}
    T = range(len(data['Time']))

    # copy over all unmapped ones
    for k in data.keys():
        mapped_data[k] = data[k]

    # 1st pass with sum mappings
    for group in mapping.keys():
        the_mapping = mapping[group]
        mapped_data[group] = {'runs': []}

        try:
            # check if is a list -> sum
            if type(the_mapping) is list:

                for run in range(params["repetitions"]):
                    this_run = []
                    runs_not_returned_by_basico = {}
                    for t in T:
                        mapping_sum = 0
                        for spe in the_mapping:
                            try:
                                mapping_sum = mapping_sum + data[spe]['runs'][run][t]
                            except KeyError:
                                mapping_sum = mapping_sum + species_not_mapped[spe]
                                try:
                                    runs_not_returned_by_basico[spe] += [species_not_mapped[spe]]
                                except KeyError:
                                    runs_not_returned_by_basico[spe] = [species_not_mapped[spe]]
                        this_run.append(mapping_sum)
                    for spe in runs_not_returned_by_basico:
                        try:
                            mapped_data[spe]['runs'].append(runs_not_returned_by_basico[spe])
                        except KeyError:
                            mapped_data[spe] = {'runs': [runs_not_returned_by_basico[spe]]}
                    mapped_data[group]['runs'].append(this_run)

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
