import sys, os
import matplotlib.pyplot as plt
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *

Chemical, Promoter, Protein = BaseSpecies(3)

Chemical + Promoter.inactive >> Promoter.active [2]
Promoter.active >> Chemical + Promoter.inactive [1]
Zero >> Protein [lambda : f'3*{Promoter.inactive}/{Promoter}']
Protein >> Zero [1]
