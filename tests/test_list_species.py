"""Tests for the List_Species container."""

from __future__ import annotations

import pytest

from mobspy import BaseSpecies
from mobspy.dsl.list_species import List_Species
from mobspy.exceptions import ValidationError


@pytest.fixture()
def species_abc():
    A, B, C = BaseSpecies(3)
    return A, B, C


class TestConstruction:
    def test_from_list(self, species_abc):
        A, B, C = species_abc
        ls = List_Species([A, B, C])
        assert len(ls) == 3

    def test_from_set(self, species_abc):
        A, B, C = species_abc
        ls = List_Species({A, B, C})
        assert len(ls) == 3

    def test_invalid_item_raises(self):
        with pytest.raises(ValidationError):
            List_Species(["not_a_species"])


class TestAppend:
    def test_append_species(self, species_abc):
        A, B, C = species_abc
        ls = List_Species([A])
        ls.append(B)
        assert len(ls) == 2

    def test_append_non_species_raises(self, species_abc):
        A, _, _ = species_abc
        ls = List_Species([A])
        with pytest.raises(ValidationError):
            ls.append("not_a_species")


class TestOperators:
    def test_or_with_species(self, species_abc):
        A, B, _ = species_abc
        ls = List_Species([A])
        ls = ls | B
        assert len(ls) == 2

    def test_or_with_list_species(self, species_abc):
        A, B, C = species_abc
        ls1 = List_Species([A])
        ls2 = List_Species([B, C])
        ls1 = ls1 | ls2
        assert len(ls1) == 3

    def test_or_invalid_raises(self, species_abc):
        A, _, _ = species_abc
        ls = List_Species([A])
        with pytest.raises(ValidationError):
            ls | "invalid"

    def test_str(self, species_abc):
        A, _, _ = species_abc
        ls = List_Species([A])
        result = str(ls)
        assert isinstance(result, str)


class TestContainerOps:
    def test_getitem(self, species_abc):
        A, B, _ = species_abc
        ls = List_Species([A, B])
        assert ls[0] is A
        assert ls[1] is B

    def test_setitem(self, species_abc):
        A, B, C = species_abc
        ls = List_Species([A, B])
        ls[1] = C
        assert ls[1] is C

    def test_is_in(self, species_abc):
        A, B, C = species_abc
        ls = List_Species([A, B])
        assert ls.is_in(A) is True
        assert ls.is_in(C) is False

    def test_insert(self, species_abc):
        A, B, C = species_abc
        ls = List_Species([A, C])
        ls.insert(1, B)
        assert len(ls) == 3
        assert ls[1] is B

    def test_extend(self, species_abc):
        A, B, C = species_abc
        ls1 = List_Species([A])
        ls2 = List_Species([B, C])
        ls1.extend(ls2)
        assert len(ls1) == 3

    def test_remove(self, species_abc):
        A, B, _ = species_abc
        ls = List_Species([A, B, A])
        removed = ls.remove(A)
        assert removed == 2
        assert len(ls) == 1

    def test_remove_repeated_elements(self, species_abc):
        A, B, _ = species_abc
        ls = List_Species([A, B, A, B])
        deduped = ls.remove_repeated_elements()
        assert len(deduped) == 2

    def test_iter(self, species_abc):
        A, B, _ = species_abc
        ls = List_Species([A, B])
        items = list(ls)
        assert len(items) == 2
