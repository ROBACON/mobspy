from mobspy import *

if __name__ == '__main__':

    """
        This is the Tree model from the paper
        We have a population of Trees
        The Trees can die, age, have different colors and be in two different forests
        The colors can change randomly from time to time
        All old Trees can reproduce, but the Three is born green and young 
    """
    Ager, Mortal, Colored, Location = BaseSpecies()
    Colored.green, Colored.yellow, Colored.brown
    Location.dense, Location.sparse
    Ager.young >> Ager.old[1 / 10 / u.year]
    Mortal >> Zero[lambda r1: 0.1/ u.year if r1.old else 0]
    Tree = Ager * Colored * Mortal * Location

    # replication
    Tree.old >> Tree + Tree.green.young[0.1/u.year]

    # competition
    Tree.dense.old + Tree.dense.young >> Tree.dense.old [1e-10 * u.decimeter**2 / u.year]

    # color cycling
    colors = ['green', 'yellow', 'brown']
    for color, next_color in zip(colors, colors[1:] + colors[:1]):
        Tree.c(color) >> Tree.c(next_color)[10/u.year]

    # initial conditions
    Tree.dense(50), Tree.dense.old(50), Tree.sparse(50), Tree.sparse.old(50)
    MySim = Simulation(Tree)
    MySim.simulation_method = 'stochastic'
    MySim.save_data = False
    MySim.plot_data = False
    MySim.duration = 100*u.years
    MySim.unit_x = 'year'
    MySim.volume = 1*u.meter**2
    MySim.repetitions = 3
    MySim.run()
    MySim.plot_stochastic(Tree.dense, Tree.sparse, Tree.brown)