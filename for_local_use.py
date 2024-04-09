from mobspy import *
import matplotlib.pyplot as plt

if __name__ == '__main__':

    def test_parameter_fit_with_units():
        A = BaseSpecies()
        A >> Zero[3 / u.hour]
        A(100)
        S1 = Simulation(A)
        S1.duration = 3 * u.hour
        S1.run(plot_data=False, level=-1)

        A = BaseSpecies()
        p = ModelParameters(1 / u.hour)
        A >> Zero[p]
        A(100)
        S2 = Simulation(A)
        S2.run(plot_data=False, level=-1)
        S2.load_experiment_data(S1.results)
        basiCO_parameter_estimation(S2, [p], verbose=False)

        print(p.value)

        assert 2.5 / u.hour <= p.value <= 3.5 / u.hour
    test_parameter_fit_with_units()
    exit()


