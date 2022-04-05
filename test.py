from mobspy import *

# TODO Plot has random order for species names

if __name__ == '__main__':

    Ager, Mortal, Colored, Location = BaseSpecies(4)
    Colored.green, Colored.yellow, Colored.brown
    Location.dense, Location.sparse
    Ager.young >> Ager.old[1 / 10 / u.year]
    Mortal >> Zero[lambda r1: 1/ u.year if r1.old else 0.01/ u.year]
    Tree = Ager * Colored * Mortal * Location

    # replication
    Tree >> Tree + Tree.green.young[0.1/u.year]

    # color cycling
    colors = ['green', 'yellow', 'brown']
    for color, next_color in zip(colors, colors[1:] + colors[:1]):
        Tree.c(color) >> Tree.c(next_color)[10/u.year]

    # competition
    Tree.dense.old + Tree.dense.young >> Tree.dense.old [0.001*u.decimeter/u.year]

    # initial conditions
    Tree.dense(50), Tree.dense.old(50), Tree.sparse(50), Tree.sparse.old(50)
    MySim = Simulation(Tree)
    MySim.simulation_method = 'stochastic'
    MySim.save_data = False
    MySim.plot_data = False
    MySim.duration = 100*u.years
    MySim.unit_x = 'year'
    MySim.compile()
    exit()
    MySim.plot_stochastic(Tree.dense, Tree.sparse, Tree.brown)


