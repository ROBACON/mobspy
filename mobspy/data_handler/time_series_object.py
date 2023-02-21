import pandas as pd
from mobspy.modules.meta_class import *


class MobsPyTimeSeries:

    def __init__(self, unprocessed_data_list):
        self.ts_data = [data['data'] for data in unprocessed_data_list]
        self.ts_parameters = [data['params'] for data in unprocessed_data_list]
        self.ts_mappings = [data['mappings'] for data in unprocessed_data_list]
        self.time_series_number = 0

    def add_ts_to_data(self, time_series):
        if isinstance(time_series, MobsPyTimeSeries):
            for ts in time_series:
                self.ts_data += ts
            for ts_parameters in self.ts_parameters:
                self.ts_parameters += ts_parameters
            for ts_mappings in self.ts_mappings:
                self.ts_mappings += ts_mappings
        elif type(time_series) == dict:
            self.ts_data += [time_series]
            self.ts_parameters += [None]
            self.ts_mappings += [None]

    def __len__(self):
        return len(self.ts_data)

    def __str__(self):
        tr = ''
        for i, data in enumerate(self.ts_data):
            tr = tr + f'Time Series {i}: \n'
            tr = tr + str(pd.DataFrame.from_dict(data)) + '\n\n'
        return tr

    def __getitem__(self, item):
        if type(item) == int:
            return self.ts_data[item]
        elif type(item) == str:
            try:
                return self.ts_data[self.time_series_number][item]
            except KeyError:
                return self._sum_reacting_species_data(item, self.time_series_number)
        elif isinstance(item, Species):
            return self.ts_data[self.time_series_number][item.get_name()]
        elif isinstance(item, Reacting_Species):
            return self._sum_reacting_species_data(item, self.time_series_number)

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
                    if reactant['characteristics'].union(reactant['object'].get_name()).issubset(set(key.split('.'))):
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

    def __getattr__(self, item):
        if 'sts' in item:
            self.time_series_number = int(item.split('_')[-1])
            return self
        elif 'runs' in item:
            return _MobsPyTimeSeriesRuns(self)
        else:
            super().__getattribute__(item)

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


class _MobsPyTimeSeriesRuns:

    def __init__(self, mobspy_ts):
        self.mobspy_ts = mobspy_ts

    def __getitem__(self, item):
        if type(item) == str:
            return [self.mobspy_ts.ts_data[i][item] for i in range(len(self.mobspy_ts))]
        elif isinstance(item, Species):
            return [self.mobspy_ts.ts_data[i][item.get_name()] for i in range(len(self.mobspy_ts))]
        elif isinstance(item, Reacting_Species):
            return [self.mobspy_ts._sum_reacting_species_data(item, i) for i in range(len(self.mobspy_ts))]


