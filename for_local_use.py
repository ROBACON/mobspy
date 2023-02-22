from mobspy import *


if __name__ == '__main__':

    B, C = BaseSpecies(2)
    C.c1, C.c2
    A = New(C)
    A.a1, A.a2
    B.b1, B.b2
    A >> Zero[2]
    B >> Zero[1]

    A(200), B(200)
    S1 = Simulation(A | B)

    S1.step_size = 0.01
    S1.save_data = False
    # S1.plot_data = False
    S1.repetitions = 3
    # MySim.initial_duration = 60
    S1.simulation_method = "deterministic"

    # MySim.duration = (A <= 0) & (B <= 0)
    S1.duration = 10

    # A(200), B(200)
    # S2 = Simulation(A | B)
    # S2.duration = 0.01

    # SC = S1 + S2
    S1.run()
    print(S1.results[A])

    # R1 = -5*A + B > 5
    # R2 = A.a1 == 0
    # R3 = B.b1 == 0

    # C1 = R1 & (R2 | R3)
    # print(C1.generate_string_from_vec_space(MySim._species_for_sbml))





