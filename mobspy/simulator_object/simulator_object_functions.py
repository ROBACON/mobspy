
def sim_remove_reaction(sim, reaction, Simulation_Constructor):
    new_sim = Simulation_Constructor(sim.model)
    new_sim._reactions_set.remove(reaction)
    return new_sim

