"""List_Species container for grouping multiple meta-species."""

from __future__ import annotations

from collections.abc import Generator, Sequence
from typing import TYPE_CHECKING, Self

from mobspy.exceptions import ValidationError

if TYPE_CHECKING:
    from mobspy.modules.species import Species


class List_Species:
    """Store species after the ``|`` operation.

    Creates a list of species that can be looped through or given
    to the simulator.

    :param list_of_species: Meta-species list to store
    """

    def __init__(self, iterable: Sequence[Species] | set[Species]) -> None:
        """Construct from an iterable of Species.

        :param iterable: Species to store
        """
        from mobspy.modules.species import Species

        self._list_species: list[Species] = []
        for item in iterable:
            if isinstance(item, Species):
                self._list_species.append(item)
            else:
                raise ValidationError("Only Species can used to construct List_Species")

    def append(self, species: Species) -> None:
        """Add species to the list.

        :param species: Meta-species to be added
        """
        from mobspy.modules.species import Species

        if not isinstance(species, Species):
            raise ValidationError("Only Species can be appended")
        else:
            self._list_species.append(species)

    def __str__(self) -> str:
        """String representation: list of meta-species names."""
        to_return = []
        for spe in self:
            to_return.append(spe.get_name())
        return str(to_return)

    def __or__(self, other: Species | List_Species) -> Self:
        """Implementation of the ``|`` operator.

        :param other: Species or List_Species to combine
        """
        from mobspy.modules.species import Species

        if isinstance(other, Species):
            self._list_species.append(other)
        elif isinstance(other, List_Species):
            self._list_species = self._list_species + other._list_species
        else:
            raise ValidationError("Operator must only be used in Species on List_Species")

        return self

    def __iter__(self) -> Generator[Species, None, None]:
        for spe in self._list_species:  # noqa: UP028
            yield spe

    def __len__(self) -> int:
        return len(self._list_species)

    def is_in(self, item: Species) -> bool:
        return item in self._list_species

    def __setitem__(self, index: int, item: Species) -> None:
        self._list_species[index] = item

    def __getitem__(self, item: int) -> Species:
        return self._list_species[item]

    def remove_repeated_elements(self) -> List_Species:
        """Create a List_Species with no repetitions.

        The order of elements is lost.

        :return: List_Species with no repeated elements
        """
        set_species = set(self._list_species)
        return List_Species(set_species)

    def insert(self, index: int, item: Species) -> None:
        self._list_species.insert(index, item)

    def extend(self, other: List_Species) -> None:
        if isinstance(other, List_Species):
            self._list_species = self._list_species + other._list_species

    def remove(self, value: object) -> int:
        """Remove all instances of the meta-species.

        :return: The number of instances removed
        """
        indexes = []
        for i, e in enumerate(self._list_species):
            if e == value or str(e) == str(value):
                indexes.append(-len(self._list_species) + i)

        for i in indexes:
            del self._list_species[i]
        return len(indexes)
