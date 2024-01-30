from mobspy import *

if __name__ == '__main__':

    def test_average_value():
        E = BaseSpecies(1)
        Zero >> E[12]
        E >> Zero[25]

        MySim = Simulation(E)
        MySim.save_data = False
        MySim.plot_data = False
        MySim.level = -1
        MySim.run()
    test_average_value()








