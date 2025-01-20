import mobspy.simulation_logging.log_scripts as simlog


def sim_remove_reaction(sim, reaction, Simulation_Constructor):
    new_sim = Simulation_Constructor(sim.model)
    new_sim._reactions_set.remove(reaction)
    return new_sim


class Simulation_Utils:

    def update_model(self, *args):

        for arg in args:

            if len(arg) != 2:
                simlog.error(' In .update_model method - '
                             'Please all parameters and species changes must be in the format: \n'
                             '(name, value)')


