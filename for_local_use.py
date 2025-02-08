from mobspy import *
from mobspy import modules

if __name__ == '__main__':

    A = BaseSpecies()
    p1, p2 = ModelParameters([1], [1 / u.hour, 2 / u.hour, 3 / u.hour])
    A >> Zero[p1 * p2]
    A(100)
    S = Simulation(A)
    S.run(duration=5 * u.hour, plot_data=False, level=-1)
    print(S.results)

    exit()

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
    Mortal >> Zero [lambda r1: 0.1 / u.year if r1.old else 0]
    Tree = Ager * Colored * Mortal * Location

    # replication
    Tree.old >> Tree + Tree.green.young[0.1/u.year]

    # competition
    Tree.dense.old + Tree.dense.young >> Tree.dense.old[1e-10 * u.decimeter ** 2 / u.year]

    # reproduction
    bf = 1e-10 * u.decimeter ** 2 / u.year
    rep_r = lambda t1, t2: 5 * bf if (Location(t1) == Location(t2) and Colored(t1) == Colored(t2)) else bf
    2 * Tree >> 2 * Tree + Tree.yound[rep_r]

    # color cycling
    colors = ['green', 'yellow', 'brown']
    for color, next_color in zip(colors, colors[1:] + colors[:1]):
        Tree.c(color) >> Tree.c(next_color)[10/u.year]

    # initial conditions
    Tree.dense(50), Tree.dense.old(50), Tree.sparse(50), Tree.sparse.old(50)
    # initial conditions
    Tree.dense(50), Tree.dense.old(50), Tree.sparse(50), Tree.sparse.old(50)
    MySim = Simulation(Tree)
    MySim.run(volume = 1 * u.meter ** 2, unit_x='year',
              duration=100 * u.years, repetitions=3, output_concentration=False,
              simulation_method='stochastic', save_data=False, plot_data=False)
    MySim.plot_stochastic(Tree.dense, Tree.sparse, Tree.brown)


