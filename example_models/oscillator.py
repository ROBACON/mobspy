import sys, os
import matplotlib.pyplot as plt
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *

# TODO look at repression and promoter reactions

Decaying = BaseSpecies(1)

TetR = Decaying*New
lcl = Decaying*New
Lacl = Decaying*New

Decaying >> Zero [0.1]

