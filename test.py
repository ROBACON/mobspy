from mobspy import *
import os


if __name__ == '__main__':


    Age, Color = BaseSpecies(2)
    Age.young, Age.old
    Color.red, Color.blue, Color.yellow
    Ecoli = Age * Color
    Test = frc.Specific_Species_Operator('Ecoli_dot_young_dot_red', Ecoli)
    print(Age(Test))
    exit()

    A, B, C = BaseSpecies(3)
    A(300) + B(200) >> C [0.1]
    Simulation(A | B | C).run()
    exit()


    def infection(r1):
        if r1.old:
            return 0.02
        else:
            return 0.01

    Age, Infection, Virus = BaseSpecies(3)
    Age.young >> Age.old [0.5]
    Default[Infection.not_infected + Virus >> Infection.infected [infection]]
    Bacteria = Age*Infection
    Bacteria(100), Virus(300)
    my_simulation = Simulation(Bacteria | Virus).compile()
    # my_simulation.parameters['simulation_method'] = 'stochastic'
    # my_simulation.parameters['plot_data'] = False
    # my_simulation.run()
    # my_simulation.plot_stochastic()
    exit()



    # def hi(**kargs):
    #   a = kargs.get('a')
    #    b = kargs.get('b')
    #    print(a, b)

    # hi(a=10, b=None)