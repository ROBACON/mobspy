import cProfile
from mobspy import *

def profile_imports():
    # Import the modules you want to profile
    import mobspy


if __name__ == "__main__":
    # cProfile.run("profile_imports()", sort='cumulative')

    A, B, C, D, E, F = BaseSpecies(6)

    A + B >> Zero[1]

    A(50), B(50), C(0)
    S = Simulation(A | B | C | D | E | F)
    S.plot_data = False
    S.level = -1

    with S.event_time(0):
        F(1)

    with S.event_condition((A <= 1) & (B <= 1)):
        C(1)

    with S.event_condition((A <= 1) & (B <= 1)):
        D(1)

    with S.event_condition(B <= 1):
        E(1)

    S.duration = 5
    print(S.compile())



