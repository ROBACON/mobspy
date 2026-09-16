"""Plot configuration must reject mistakes and own mutable style dictionaries."""

import pytest

from mobspy import BaseSpecies, Simulation
from mobspy.exceptions import ValidationError


@pytest.mark.parametrize("composed", [False, True])
def test_plot_styles_are_copied(composed):
    A = BaseSpecies()
    first = Simulation(A)
    target = first + Simulation(A) if composed else first
    style = {"A": {"color": "red", "label": "Population"}}

    target.add_plot_params(style, title="Growth")
    style["A"]["color"] = "blue"

    assert target.plot_parameters["A"]["color"] == "red"
    assert target.plot_parameters["title"] == "Growth"
    if composed:
        assert "A" not in first.plot_parameters


@pytest.mark.parametrize("composed", [False, True])
def test_invalid_plot_argument_fails_without_partial_update(composed):
    A = BaseSpecies()
    first = Simulation(A)
    target = first + Simulation(A) if composed else first
    original = dict(target.plot_parameters)

    with pytest.raises(ValidationError, match="expects dictionaries"):
        target.add_plot_params({"title": "Changed"}, A, color="red")

    assert target.plot_parameters == original
