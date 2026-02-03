"""Simple tests for exponents in MobsPy."""
from mobspy import BaseSpecies, Simulation, u, New

infection_by_antibio_sender = 3e-11 *u.mL / u.min
k_tet = 1208585 *u.counts / u.mL
n_tet = 1.0
qty_sender = 1000 *u.counts / u.mL
qty_tet = 1e10 *u.counts / u.mL

def test():

    def antibio_uptake_senders(tet, b):
        #breakpoint()
        mu_monod = (infection_by_antibio_sender)*(k_tet**n_tet/((k_tet**n_tet)+qty_sender**n_tet))
        rate_monod = mu_monod * b
        #test_rate = (infection_by_antibio_sender)*(k_tet/((k_tet)+qty_sender))
        return rate_monod #test_rate

    # Base species definition
    Antibiotic = BaseSpecies()

    Tet = New(Antibiotic)
    Tet(qty_tet)

    tet = BaseSpecies()
    tet.tetF, tet.tetT

    Bacterium = tet
    Sender = New(Bacterium)

    Sender(qty_sender)

    # Reaction
    Tet + Sender.tetF >> Sender.tetT [antibio_uptake_senders] 

    return Simulation( Sender | Tet )


if __name__ == "__main__":
    S = test()
    S.run(duration=1222 * u.min, unit_x=u.min, unit_y=u.counts / u.mL, plot_data=False)
