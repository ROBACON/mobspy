from mobspy import *
from mobspy.modules.ode_operator import dt
from testutils import compare_model, compare_model_ignore_order


if __name__ == "__main__":

    def ode_neg_test():
        Neg, NegR = BaseSpecies()
        NegR.comp1

        dt[Neg] >> -Neg
        dt[NegR] >> -NegR.comp1

        S = Simulation(Neg | NegR)
        assert compare_model(S.compile(), "test_tools/model_ode_neg_test.txt")
    ode_neg_test()

    exit()

    def ode_with_meta_species_multiplication():
        """
        Test ODE syntax with meta-species multiplication (Cartesian product).

        Models a system where Location and State are combined, and decay rate
        depends on the specific combination via lambda rate function.

        Species generated: A.here.alive, A.here.dead, A.there.alive, A.there.dead
        Only alive species decay, with location-dependent rates.
        """
        Location, State = BaseSpecies(2)
        Location.here, Location.there
        State.alive, State.dead

        A = Location * State

        # ODE: dA/dt = -k*A where k depends on location
        # here decays faster than there
        dt[A.alive] >> -A.alive  [lambda r1: 0.2 if Location.here else 0.1]

        A.here.alive(100), A.there.alive(100)
        S = Simulation(A)
        print(S.compile())
    # ode_with_meta_species_multiplication()
