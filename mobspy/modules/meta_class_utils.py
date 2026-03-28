"""This model stores function used by the meta_class.py module"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from mobspy.exceptions import ValidationError


def count_stoichiometry(
    lst: Sequence[tuple[float, str] | str],
) -> dict[str, float | int]:
    to_return: dict[str, float | int] = {}

    for e in lst:
        stoichiometry, species = (1, e) if isinstance(e, str) else e
        if species not in to_return:
            to_return[species] = 0

        to_return[species] += stoichiometry

    return to_return


def combine_references(species1: Any, species2: Any) -> set[Any]:
    """Combine the sets of references of two species

    :param species1: (Meta-species object)
    :param species2: (Meta-species object)
    """
    return species1.get_references().union(species2.get_references())


def check_orthogonality_between_references(references: set[Any]) -> None:
    """Check if meta-species objects inside a reference do
    not have characteristics in common. The sets of
    characteristics directly added to species must be
    independent.

    :param references: (set) set of meta-species objects
        to check for independence
    :raise simlog.error: raises error if there are repeated
        characteristics in different meta-species
    """
    for i, reference1 in enumerate(references):
        for j, reference2 in enumerate(references):
            if i == j:
                continue

            if (
                len(
                    reference1.get_characteristics().intersection(
                        reference2.get_characteristics()
                    )
                )
                != 0
            ):
                raise ValidationError(
                    "The same characteristic can only be "
                    "shared through inheritance. " + "There are two characteristics "
                    "directly added to two "
                    "meta-species \n"
                    f"Repetition in: "
                    f"{reference1}, {reference2}"
                    f"Characteristics: "
                    f"{reference1.get_characteristics()}"
                    f", "
                    f"{reference2.get_characteristics()}"
                )


def unite_characteristics(species: list[Any] | None) -> set[str]:
    """This function unites the characteristics of all the given species

    :param species: (list of species or List_Species object)
    """
    characteristics: set[str] = set()

    if species is not None:
        for spe in species:
            characteristics = characteristics.union(spe.get_characteristics())

    return characteristics


def create_orthogonal_vector_structure(
    species: list[Any] | set[Any],
) -> dict[str, Any]:
    """This creates the independent state-structure for the
    model. It is just a dictionary where the keys are
    characteristics and the values are meta-species objects
    that have been directly added to that object (no
    inheritance or New).
    It simplifies the code by allowing to easily keep track
    of the 'axis' of each characteristic. Allowing for easy
    transformation on the products and others.

    :param species: (meta-species objects) - meta-species
        objects used in a model

    :return ref_characteristics_to_object: (dict) a
        dictionary where the keys are characteristics and
        the values are meta-species objects that have been
        directly added to that object
    """
    ref_characteristics_to_object: dict[str, Any] = {}
    for spe in species:
        for prop in spe.get_references():
            for cha in prop.get_characteristics():
                if cha not in ref_characteristics_to_object:
                    ref_characteristics_to_object[cha] = prop
                elif ref_characteristics_to_object[cha] == prop:
                    pass
                else:
                    raise ValidationError(
                        "The same characteristic can only "
                        "be shared through inheritance. "
                        + "There are two characteristics "
                        "directly added to two "
                        "meta-species \n"
                        f"Repetition in: {spe}, {ref_characteristics_to_object[cha]} \n"
                        f"Characteristics: {spe.get_characteristics()}, "
                        f"{ref_characteristics_to_object[cha].get_characteristics()} \n"
                    )

    return ref_characteristics_to_object


if __name__ == "__main__":
    pass
