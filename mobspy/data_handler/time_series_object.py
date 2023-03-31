import pandas as pd
from mobspy.modules.meta_class import *


class MobsPyTimeSeries:

    def __init__(self, data_dict):
        self.ts_data = data_dict['data']
        self.ts_parameters = data_dict['params']
        self.ts_models = data_dict['models']

    def to_dict(self):
        return {'data': self.ts_data,
                'params': self.ts_parameters,
                'models': self.ts_models}

    def __len__(self):
        return len(self.ts_data)

    def __str__(self):
        tr = ''
        for i, data in enumerate(self.ts_data):
            tr = tr + f'Time Series {i}: \n'
            tr = tr + str(pd.DataFrame.from_dict(data)) + '\n\n'
        return tr

    def add_ts_to_data(self, time_series):
        if type(time_series) == dict:
            self.ts_data += [time_series]

    def __getitem__(self, item):
        to_return = []
        if type(item) == int:
            return self.ts_data[item]
        elif type(item) == str:
            try:
                for ts in self:
                    to_return.append(ts[item])
            except KeyError:
                for i in range(len(self)):
                    to_return.append(self._sum_reacting_species_data(item, i))
        elif isinstance(item, Species):
            for ts in self:
                to_return.append(ts[item.get_name()])
        elif isinstance(item, Reacting_Species):
            for i in range(len(self)):
                to_return.append(self._sum_reacting_species_data(item, i))

        if len(to_return) == 1:
            return to_return[0]
        else:
            return to_return

    def _sum_reacting_species_data(self, item, ts_index):

        def _sum_element_by_element(l1, l2):
            rt = []
            for e1, e2 in zip(l1, l2):
                rt.append(e1 + e2)
            return rt

        time_series = self.ts_data[ts_index]
        to_return = [0 for _ in range(len(time_series['Time']))]

        found_flag = False
        if isinstance(item, Reacting_Species):
            for reactant in item.list_of_reactants:
                for key in time_series:
                    reactant_name = reactant['object'].get_name()
                    reactant_full_name = reactant['characteristics'].union([reactant_name])
                    if reactant_full_name.issubset(set(key.split('.'))):
                        to_return = _sum_element_by_element(to_return, time_series[key])
                        found_flag = True
        elif type(item) == str:
            for key in time_series:
                if set(item.split('.')).issubset(set(key.split('.'))):
                    to_return = _sum_element_by_element(to_return, time_series[key])
                    found_flag = True

        if not found_flag:
            simlog.error(f'{item} was not found in data')

        return to_return

    def __iter__(self):
        for ts in self.ts_data:
            yield ts

    def get_max_time_for_species(self, species):
        if type(species) != str:
            species = species.get_name()
        max_length = 0
        max_ts = None
        for ts in self:
            if species in ts.keys():
                if len(ts) > max_length:
                    max_length = len(ts)
                    max_ts = ts
        return max_ts['Time']
