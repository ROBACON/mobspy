from mobspy import *

if __name__ == '__main__':

    A = BaseSpecies()

    A >> Zero [1]

    A(100)
    S = Simulation(A)
    S.duration = 10
    S.run()
    print(S.to_dataframe()[0])
    exit()


    #A = BaseSpecies()

    #Set[A >> 2*A [1]].at(5)
   # Set[A >> 2*A [2]].when(A > 5)

    #S = Simulation(A)
    #with S.event_time(5):
    #    A(10)

    print(S.compile())


