from mobspy import *
import numpy as np

if __name__ == '__main__':

   A = BaseSpecies()

   A >> Zero[1]
   S = Simulation(A)
   S.run()
