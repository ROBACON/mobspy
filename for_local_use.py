from mobspy import *


if __name__ == '__main__':

    def text_even_more_complex_assignments():
        Hi = BaseSpecies()
        Hi.h1, Hi.h2
        A, B, C, D = New(Hi)

        A.assign(B * ((C + 5) / D))

        S = Simulation(A | B | C | D)
        print(S.compile())
    text_even_more_complex_assignments()
    exit()


    # Fixed constant and their units
    Ein = 2.92  # ÂµE.s^-1
    K = 5.1e-9  # mL.cm^-1.cell^-1
    V = 200  # mL
    D = 3.7  # cm

    Microalga, Nutrient, Light, Complex, Space = BaseSpecies(5)
    Complex.activated, Complex.inactivated

    # Assignment
    Light.assign((1 - 10**(-K * Microalga / (V * D))) * Ein / Microalga)

    # Nutrient reaction : 1st equilibrium (units are not set)d
    Microalga + Nutrient >> Complex.inactivated[1/ u.hour]
    Complex.inactivated >> Microalga + Nutrient[1/ u.hour]

    # Light reaction : 2nd equilibrium (units are not set)
    Complex.inactivated + Light >> Complex.activated[1/u.hour]
    Complex.activated >> Complex.inactivated + Light[1/u.hour]

    # Mitosis (units are not set)
    Complex.activated + Space >> Microalga + Microalga[1/ u.hour]

    # Simulation
    Microalga(0.026 * 30.1e6 * 200), Nutrient(8.7e12), Complex(0), Space(9.99e10)
    S = Simulation(Microalga | Nutrient | Complex | Light | Space)
    S.save_data = False
    S.method = 'deterministic'
    S.duration = 400 * u.hour
    S.unit_x = u.hour
    S.plot_data = False
    S.level = 1
    S.run()
    exit()

    # TODO: Fix All in Assignments

    A, B, P = BaseSpecies()

    B >> 2*B [lambda r: (100-r)*r]
    B >> Zero [1]

    A.assign(B**2/(B**2 + 50))

    A >> A + P [1]
    P >> Zero [1]

    S = Simulation(A | B | P)
    S.run()
    print(S.compile())


