"""Test in order to check the exponent in Mobspy"""
from test_exponent_model import Model, InitialQuantities
from mobspy import u

def test_exponent():
    qty_sender = 1000 * u.counts / u.mL
    qty_tet = 1e10 * u.counts / u.mL

    qty_initial = InitialQuantities(
        qty_sender=qty_sender,
        qty_tet=qty_tet,
    )
    model = Model()

    S = model.get_simulation_exp(qty_initial)
    S.run(duration=1222 * u.min, unit_x=u.min, unit_y=u.counts / u.mL, plot_data=False)

def test_default():
    qty_sender = 1000 * u.counts / u.mL
    qty_tet = 1e10 * u.counts / u.mL

    qty_initial = InitialQuantities(
        qty_sender=qty_sender,
        qty_tet=qty_tet,
    )
    model = Model()

    S = model.get_simulation_default(qty_initial)
    S.run(duration=1222 * u.min, unit_x=u.min, unit_y=u.counts / u.mL, plot_data=False)

if __name__ == "__main__":
    test_exponent()
    test_default()