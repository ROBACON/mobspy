from mobspy import *
from mobspy import modules

if __name__ == '__main__':

    Spe_Dict = {'E_1':None, 'Eb_1':None, 'D_2': None}

    L = BaseSpecies(list(Spe_Dict.keys()))

    for spe in L:
        Spe_Dict[spe.get_name()] = spe

    print(Spe_Dict)

