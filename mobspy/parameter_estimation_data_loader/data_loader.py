import mobspy.simulation_logging.log_scripts as simlog
import mobspy.data_handler.time_series_object as tso

class Experimental_Data_Holder:

    def __init__(self):
        self.experimental_data = None

    def load_experiment_data(self, data):

        flag_jump_checks = True if isinstance(data, tso.MobsPyList_of_TS) else False

        if type(data) != list and not flag_jump_checks:
            simlog.error('Data added must be in the format of list with each element being a dictionary '
                         'with species names and time as keys or a MobsPy results object')
        for e in data:
            if type(e) != dict and not flag_jump_checks:
                simlog.error('Data added must be in the format of list with each element being a dictionary '
                             'with species names and time as keys or a MobsPy results object')

        self.experimental_data = data

