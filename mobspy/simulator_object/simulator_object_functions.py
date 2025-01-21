from scipy.constants import value

import mobspy.simulation_logging.log_scripts as simlog
from mobspy.modules.mobspy_parameters import _Internal_Parameter_Constructor


def sim_remove_reaction(sim, reaction, Simulation_Constructor):
    new_sim = Simulation_Constructor(sim.model)
    new_sim._reactions_set.remove(reaction)
    return new_sim


class Simulation_Utils:

    def update_model(self, *args):

        if not self._list_of_models:
            simlog.error('In .update_model method - \n'
                         'The model was not compiled yet. The update_model method is reserved for simulations that '
                         'have already been compiled')

        for arg in args:

            if len(arg) != 2:
                simlog.error('In .update_model method - \n'
                             'Please all parameters and species changes must be in the format: \n'
                             '(name, value)')

            self._update_from_compiler(arg)

        # This does not look to be necessary
        # We need to regenerate the sbml files - each parameter has a single sbml file
        # They are generated in _assemble_multi_simulation_structure - which is a partial compilation
        # self._assemble_multi_simulation_structure()


    def _update_from_compiler(self, arg):

        if isinstance(arg[0], _Internal_Parameter_Constructor):
            self._update_parameter(arg)

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

        if type(arg) == str:
            parameter_str = arg[0]
        else:
            parameter_str = arg[0].get_name()

        # Update values on standard compiler
        for model in self._list_of_models:
            try:
                model['parameters_for_sbml'][parameter_str ] = (value_to_update, 'dimensionless')
            except KeyError:
                simlog.error(f'The parameter named {parameter_str } was not found in the model')

        # Update parameter in model_parameters in the simulation object
        try:
            if not iterable:
                self.model_parameters[parameter_str]['values'] = [arg[1]]
            else:
                self.model_parameters[parameter_str]['values'] = arg[1]
        except KeyError:
            simlog.error(f'The parameter named {parameter_str} was not found in the model')

        self.model_parameters[parameter_str]['object'].update_value(arg[1])

