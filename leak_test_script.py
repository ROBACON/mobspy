import tracemalloc
import gc
import time
from mobspy import *
from tqdm import tqdm
import matplotlib.pyplot as plt
import os
import basico
import pint
from numpy import issubdtype, floating, integer, array


if __name__ == '__main__':

    import tracemalloc
    import time

    # Start tracing memory allocations
    tracemalloc.start()


    def run_simulation():
        A = BaseSpecies()
        A >> Zero[1]
        A(100)
        S = Simulation(A)
        S.plot_data = False
        S.level = -1
        S.run()


    # Monitor memory usage over multiple iterations
    for i in range(100):  # Adjust the number of iterations as needed
        run_simulation()
        gc.collect()
        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('traceback')

        if i % 10 == 0:  # Take a dump every 10 iterations
            snapshot.dump(f"snapshot_{i}.dat")

        print(f"Iteration {i + 1} Memory Usage:")
        for stat in top_stats[:10]:
            print(stat)

    tracemalloc.stop()



