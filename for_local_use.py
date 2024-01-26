from mobspy import *

if __name__ == '__main__':

    def no_ab_model(duration=40, rate=1, init_res=10000, init_bact=1000, init_atp=0):
        Res, Bact, ATP = BaseSpecies()

        Res(init_res / u.ul)
        Bact(init_bact / u.ul)
        ATP(init_atp / u.ul)

        Res + Bact >> Bact + Bact + ATP [rate*u.ul/u.hours]

        S = Simulation(Res | Bact | ATP)
        S.compile()

        try:
            Res, Bact, ATP = BaseSpecies()

            Res(init_res / u.ul)
            Bact(init_bact / u.ul)
            ATP(init_atp / u.ul)

            Res + Bact >> Bact + Bact + ATP [rate*u.ul/u.meters]

            S = Simulation(Res | Bact | ATP)
            S.compile()
            assert False
        except:
            pass

        Res, Bact, ATP = BaseSpecies()

        Res(init_res / u.ul)
        Bact(init_bact / u.ul)
        ATP(init_atp / u.ul)

        Res + Bact >> Bact + Bact + ATP [rate*u.ul/u.hours]

        S = Simulation(Res | Bact | ATP)
        S.compile()
        assert True

    no_ab_model()





