import itertools
from copy import deepcopy


def assign_values_to_model(parameter_name, parameter_value, models, locations):

    for location in locations:
        if location == '$sbml':
            for model in models:
                try:
                    model['parameters_for_sbml'][parameter_name] = (parameter_value, 'dimensionless')
                except KeyError:
                    pass
        else:
            for model in models:
                try:
                    model['species_for_sbml'][location] = parameter_value
                    model['species_not_mapped'][location.replace('_dot_', '.')] = parameter_value
                except KeyError:
                    pass


def generate_all_sbml_models(model_parameters, list_of_models):

    names = []
    used_in = []
    values = []

    to_return = []

    if model_parameters == {}:
        return [list_of_models], []

    keys = sorted(model_parameters.keys())
    for key in keys:
        item = model_parameters[key]

        names.append(item['name'])
        used_in.append(item['used_in'])
        try:
            if len(item['values']):
                values.append(item['values'])
        except TypeError:
            values.append([item['values']])

    parameter_list_of_dic = []
    for v in itertools.product(*values):
        parameter_dic = {}
        for i, name in enumerate(names):
            assign_values_to_model(name, v[i], list_of_models, used_in[i])
            parameter_dic[name] = v[i]

        parameter_list_of_dic.append(parameter_dic)
        to_return.append(deepcopy(list_of_models))

    return to_return, parameter_list_of_dic


def unite_parameter_dictionaries(dict_1, dict_2):

    for key in dict_2:
        if key not in dict_1:
            dict_1[key] = dict_2[key]
        else:
            new_used_in = dict_1[key]['used_in'].union(dict_2[key]['used_in'])
            dict_1[key]['used_in'] = new_used_in

    return dict_1

