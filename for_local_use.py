from mobspy import *

if __name__ == '__main__':

    def second_test_model():

        A, B = BaseSpecies()
        A >> Zero[1]
        B >> Zero[2]
        A(100), B(200)
        S1 = Simulation(A)
        S1.run(plot_data=False)

        A, B = BaseSpecies()
        A >> Zero[1]
        B >> Zero[2]
        A(100), B(200)
        S2 = Simulation(A)
        S2.run(plot_data=False)
        exp = [S1.results.return_pandas()[0], S2.results.return_pandas()[0]]

        A, B = BaseSpecies()
        a = ModelParameters(0.01)
        A >> Zero[a]
        B >> Zero[2*a]
        A(100)
        S3 = Simulation(A)
        results = basiCO_parameter_estimation(S3, experimental_data=exp, parameters_to_estimate=[a])
        return results
    second_test_model()

    def first_test_model():

        A = BaseSpecies()
        A >> Zero [1]
        A(100)
        S1 = Simulation(A)
        S1.run(level=-1, plot_data=False)

        A = BaseSpecies()
        k = ModelParameters(0.5)

        A >> Zero[k]
        A(100)
        S2 = Simulation(A)
        S2.load_experiment_data(S1.results)
        results = basiCO_parameter_estimation(S2)
        return results









