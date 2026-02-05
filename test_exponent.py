"""Simple tests for exponents in MobsPy."""

from mobspy import BaseSpecies, Simulation, u

default_rate = 3e-11 / u.min
k_coeff = 1208585 * u.counts / u.mL
n_hill = 1.0
init_a_qty = 1000 * u.counts / u.mL


def test():
    def antibio_uptake_senders(a):
        mu_monod = (default_rate) * (
            k_coeff**n_hill / ((k_coeff**n_hill) + init_a_qty**n_hill)
        )
        rate_monod = mu_monod * a
        return rate_monod

    # Base species definition
    A = BaseSpecies()
    A(init_a_qty)

    # Reaction
    A >> 2 * A[antibio_uptake_senders]

    S = Simulation(A)
    S.duration = 1222 * u.min
    S.unit_x = u.min
    S.unit_y = u.counts / u.mL
    S.plot_data = False
    S.save_data = False
    S.run()


def test_rate_function_with_exp():
    def get_model():
        default_rate_ = 3e-11 / u.min
        qty_A = 1000 * u.counts / u.mL
        K = 1e10 * u.counts / u.mL
        n = 1.0

        def antibio_uptake_senders(a):
            rate_cst = default_rate_ * (K**n / (K**n + qty_A**n))
            rate = rate_cst * a
            return rate

        A = BaseSpecies()
        A(qty_A)

        A >> 2 * A[antibio_uptake_senders]
        return Simulation(A)

    S = get_model()
    S.run(duration=1222 * u.min, unit_x=u.min, unit_y=u.counts / u.mL, plot_data=False)


if __name__ == "__main__":
    test()
    test_rate_function_with_exp()
