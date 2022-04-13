from mobspy import *

# TODO Plot has random order for species names

if __name__ == '__main__':

    A, B, C, D = BaseSpecies(4)
    for i in range(10):
        A.c(f's_{i}') + B >> A.c(f's_{i + 1}') [1]
        if i == 5:
            A.c(f's_{i}') + C >> Zero [1]
        if i > 5:
            A.c(f's_{i}') + D >> 2*A [1]

    exit()

    BindingSite = BaseSpecies(1)
    BindingSite.n_0, BindingSite.n_1
    A, B = New(BindingSite, 2)
    AB = A * B
    ABC = A * B * C
    Rev[A.n_0 + B.n_0 >> AB.n_1][1,1]
    A.n_1 + C.pho_0 >> ABC.n_1
    # ABC >> AB + C.pho_1
    # A + C.pho_1 >> AC
    # AC >> A + C.pho_2
    Simulation(A | B | AB).compile()

    for i in [1, 3, 5]:
        A.c(f'p_{i}') + B >> AB
    exit()

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


