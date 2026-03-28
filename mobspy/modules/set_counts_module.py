from __future__ import annotations

import sys
from typing import Any

from numpy import floating as np_float_  # noqa: E402
from numpy import integer as np_int_  # noqa: E402
from pint import Quantity  # noqa: E402

from mobspy.exceptions import ValidationError
from mobspy.modules.meta_class import (  # noqa: E402
    List_Species,
    Reacting_Species,
    Species,
)
from mobspy.modules.mobspy_parameters import (  # noqa: E402
    Internal_Parameter_Constructor as mp_Mobspy_Parameter,
)


def set_counts(count_dic: dict[Any, Any]) -> List_Species:
    """Adds counts to meta-species using a given dictionary.

    Keys from this dictionary can be either meta-species objects or
    strings. Items must be the assigned counts to species.

    :param count_dic: Dictionary where keys can be either
        meta-species or strings
    :raise simlog.error:
        - If a count is assigned to a Reacting_Species with more
          than one meta-species.
        - If there are two species with the same name and a string
          assignment is performed.
        - If the keys are not strings or meta-species.
        - If the counts are not Quantities, Floats or Ints.
        - If the species was not found in the stack.
    :return: List_Species object. All meta-species that had a count
        assigned in this dictionary will be returned as a
        List_Species which can be passed as a model to the
        simulation object
    """
    new_count_dict: dict[Any, Any] = {}
    for key, item in count_dic.items():
        if (
            type(item) == int  # noqa: SIM101, E721
            or type(item) == float  # noqa: E721
            or isinstance(item, Quantity)
            or isinstance(item, mp_Mobspy_Parameter)
        ):
            new_count_dict[key] = item
        elif isinstance(item, (np_int_, np_float_)):
            new_count_dict[key] = float(item)
        else:
            raise ValidationError(
                "Reactant_species count assignment does not"
                f" support the type {type(item)}"
            )
    count_dic = new_count_dict

    def find_species() -> set[Species]:
        found_species: set[Species] = set()

        frame = sys._getframe(1)
        while frame is not None:
            for ns in (frame.f_locals, frame.f_globals):
                for _name, obj in ns.items():
                    try:
                        if isinstance(obj, Species) and type(obj) != type:  # noqa: E721
                            found_species.add(obj)
                    except AttributeError:
                        pass
            frame = frame.f_back

        return found_species

    all_found_species: set[Species] = set()
    for key in count_dic:
        if type(key) == str:  # noqa: E721
            all_found_species = find_species()

    model: set[Any] = set()
    for key, item in count_dic.items():
        if type(key) == str:  # noqa: E721
            already_found = False
            str_name = key.split(".")[0]
            str_characteristics = set(key.split(".")[1:])
            for spe in all_found_species:
                if spe.get_name() == str_name and not already_found:
                    already_found = True
                    temp_set = set(str_characteristics)
                    temp_set.discard("all$")
                    if temp_set.issubset(spe.get_all_characteristics()):
                        spe.add_quantities(str_characteristics, item)
                    else:
                        raise ValidationError(
                            "Characteristics not found in species with equal name"
                        )
                    model.add(spe)
                elif spe.get_name() == str_name and already_found:
                    raise ValidationError(
                        "There are two different meta-species with"
                        " the same name. Set_counts cannot resolve"
                    )
            if not already_found:
                raise ValidationError(
                    f"Meta-species with the following name {key} not found"
                )
        else:
            try:
                if isinstance(key, Species) or isinstance(key, Reacting_Species):  # noqa: SIM101
                    if not isinstance(key, Species):
                        if len(key.list_of_reactants) != 1:
                            raise ValidationError(
                                "Assignment used incorrectly."
                                " Only one species at a time"
                            )
                        model.add(key.list_of_reactants[0]["object"])
                    if isinstance(key, Species):
                        model.add(key)
                    key(item)
            except AttributeError as e:
                raise ValidationError(
                    "Keys must be either meta-species or strings"
                ) from e

    return List_Species(model)
