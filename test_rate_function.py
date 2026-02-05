"""Test in order to check the exponent in Mobspy"""

# from test_exponent_model import Model, InitialQuantities
from mobspy import u, BaseSpecies, Simulation


def test_rate_function_right_multiplication():
    def get_model():
        default_rate = 3e-11 / u.min

        def antibio_uptake_senders(a):
            rate = default_rate * a
            return rate

        A = BaseSpecies()
        A(1000 * u.counts / u.mL)

        A >> 2 * A[antibio_uptake_senders]
        return Simulation(A)

    S = get_model()
    S.run(duration=1222 * u.min, unit_x=u.min, unit_y=u.counts / u.mL, plot_data=False)


def test_rate_function_inside_multiplication():
    def get_model():
        def antibio_uptake_senders(a):
            rate = (3e-11 / u.min) * a
            return rate

        A = BaseSpecies()
        A(1000 * u.counts / u.mL)

        A >> 2 * A[antibio_uptake_senders]
        return Simulation(A)

    S = get_model()
    S.run(duration=1222 * u.min, unit_x=u.min, unit_y=u.counts / u.mL, plot_data=False)


def test_rate_function_left_multiplication():
    def get_model():
        default_rate = 3e-11 / u.min

        def antibio_uptake_senders(a):
            rate = a * default_rate
            return rate

        A = BaseSpecies()
        A(1000 * u.counts / u.mL)

        A >> 2 * A[antibio_uptake_senders]
        return Simulation(A)

    S = get_model()
    S.run(duration=1222 * u.min, unit_x=u.min, unit_y=u.counts / u.mL, plot_data=False)


if __name__ == "__main__":
    test_rate_function_right_multiplication()
    test_rate_function_inside_multiplication()
    test_rate_function_left_multiplication()
