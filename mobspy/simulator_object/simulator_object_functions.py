from pint.registry import Quantity
from scipy.constants import value

import mobspy.simulation_logging.log_scripts as simlog
from mobspy.modules.mobspy_parameters import Internal_Parameter_Constructor
from mobspy.modules.species_string_generator import construct_species_char_list as sp_construct_species_char_list
from mobspy.modules.species_string_generator import construct_all_combinations as sp_construct_all_combinations
from mobspy.modules.unit_handler import convert_counts as uh_convert_counts


def sim_remove_reaction(sim, reaction, Simulation_Constructor):
    new_sim = Simulation_Constructor(sim.model)
    new_sim._reactions_set.remove(reaction)
    return new_sim


class Simulation_Utils:

    def update_model(self, *args):

        # Check if the model was already compiled
        if not self._list_of_models:
            simlog.error('In .update_model method - \n'
                         'The model was not compiled yet. The update_model method is reserved for simulations that '
                         'have already been compiled')

        # Check every argument to see if it is in the proper format - len 2
        for arg in args:

            if len(arg) != 2:
                simlog.error('In .update_model method - \n'
                             'Please all parameters and species changes must be in the format: \n'
                             '(name, value)')

            self._update_from_compiler(arg)


    def _update_from_compiler(self, arg):

        try:
            is_species = arg[0].is_spe_or_reac()
        except:
            is_species = False

        if isinstance(arg[0], Internal_Parameter_Constructor):
            self._update_parameter(arg)
        elif type(arg[0]) == str:
            test_model = self._list_of_models[0]

            # Check to see if string is in parameters
            try:
                test_model['parameters_for_sbml'][arg[0]]
                self._update_parameter(arg)
                not_parameter = False
            except KeyError:
                not_parameter = True

            # Check to see if string is in species
            try:
                test_model['species_for_sbml'][arg[0].replace('_dot_', '.')]
                self._update_species(arg)
                not_species = False
            except KeyError:
                not_species = True

            # If it is in neither - throw an error
            if not_species and not_parameter:
                simlog.error(f"The string {arg[0]} was not found either in parameters or species")

        elif is_species:
            self._update_species(arg)

        else:
            simlog.error('Placeholder error for now')


    def _update_parameter(self, arg):

        try:
            iterable = iter(arg[1])
        except TypeError:
            iterable = False

        if iterable:
            value_to_update = arg[1][0]
        else:
            value_to_update = arg[1]

        if type(arg[0]) == str:
            parameter_str = arg[0]
        else:
            parameter_str = arg[0].get_name()

        # Update values on standard compiler
        for model in self._list_of_models:
            try:
                model['parameters_for_sbml'][parameter_str] = (value_to_update, 'dimensionless')
            except KeyError:
                simlog.error(f'The parameter named {parameter_str } was not found in the model')

        parameter_object = self.model_parameters[parameter_str]['object']
        parameter_object.update_value(arg[1])

        # Update parameter in model_parameters in the simulation object
        try:
            if not iterable:
                self.model_parameters[parameter_str]['values'] = [parameter_object.value]
            else:
                self.model_parameters[parameter_str]['values'] = parameter_object.value
        except KeyError:
            simlog.error(f'The parameter named {parameter_str} was not found in the model')


    def _update_species(self, arg):

        # Prepare count
        if 'volume' not in self.__dict__:
            volume = 1
        else:
            volume = self.__dict__['volume']
        dimension = self.__dict__['dimension']

        spe_count = uh_convert_counts(arg[1], volume, dimension)

        # Get query - construct all combinations - or just one
        query = arg[0].get_query_characteristics()
        if 'all$' in query:

            spe_string_list = sp_construct_all_combinations(arg[0], query,
                                                            self.orthogonal_vector_structure, symbol='_dot_')

            for spe_string in spe_string_list:
                self._list_of_models[0]['species_for_sbml'][spe_string] = spe_count

        else:

            spe_string = sp_construct_species_char_list(arg[0], query,
                                                        self.orthogonal_vector_structure, symbol='_dot_')

            self._list_of_models[0]['species_for_sbml'][spe_string] = spe_count




