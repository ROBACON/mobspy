from mobspy import *

if __name__ == '__main__':

    def test_simple_fit():

        A = BaseSpecies()
        A >> Zero [1]
        A(100)
        S1 = Simulation(A)
        S1.run(level=-1, plot_data=False, step_size=1)

        A = BaseSpecies()
        k = ModelParameters(0.5)

        A >> Zero [k]
        A(100)
        S2 = Simulation(A)
        S2.load_experiment_data(S1.results)
        basiCO_parameter_estimation(S2, [k], bound=(0, 2))

        assert 0.8 <= k.value <= 1.2
    test_simple_fit()
    exit()








