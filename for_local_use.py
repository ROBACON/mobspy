from mobspy import *


if __name__ == '__main__':

    A, B = BaseSpecies(2)
    A.a1, A.a2
    B.b1, B.b2
    MySim = Simulation(A | B)
    A >> Zero [2]
    B >> Zero [1]
    A(200)
    with MySim.event_delay(0) as _:
        A >= 0
        if A <= 0:
            A(200)
    MySim.step_size = 0.01
    MySim.save_data = False
    print(MySim.compile())
    MySim.run()

    # R1 = -5*A + B > 5
    # R2 = A.a1 == 0
    # R3 = B.b1 == 0

    # C1 = R1 & (R2 | R3)
    # print(C1.generate_string_from_vec_space(MySim._species_for_sbml))





