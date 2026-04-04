"""Tests for the mobspy.types dataclass backward compatibility."""

from __future__ import annotations

from mobspy.types import (
    AssignmentData,
    CompiledModel,
    CompilerResult,
    EventData,
    ReactionData,
    SimulationEventData,
)


class TestReactionData:
    def test_attribute_access(self):
        r = ReactionData(
            reactants=[(1, "A")],
            products=[(1, "B")],
            kinetics="k * A",
        )
        assert r.reactants == [(1, "A")]
        assert r.products == [(1, "B")]
        assert r.kinetics == "k * A"

    def test_defaults(self):
        r = ReactionData()
        assert r.reactants == []
        assert r.products == []
        assert r.kinetics == ""


class TestEventData:
    def test_attribute_access(self):
        e = EventData(
            trigger="time > 10",
            delay=5,
            assignments=[("A", 100)],
        )
        assert e.trigger == "time > 10"
        assert e.delay == 5
        assert e.assignments == [("A", 100)]


class TestAssignmentData:
    def test_attribute_access(self):
        a = AssignmentData(species="A", expression="k * B")
        assert a.species == "A"
        assert a.expression == "k * B"


class TestSimulationEventData:
    def test_attribute_access(self):
        s = SimulationEventData(
            event_time=10.0,
            event_counts=[{"species": "A", "quantity": 5}],
            trigger="true",
        )
        assert s.event_time == 10.0
        assert len(s.event_counts) == 1
        assert s.trigger == "true"


class TestCompiledModel:
    def test_dict_style_access(self):
        m = CompiledModel(
            species_for_sbml={"A": 10, "B": 0},
            reactions_for_sbml={
                "r1": ReactionData(
                    reactants=[(1, "A")],
                    products=[(1, "B")],
                    kinetics="k * A",
                )
            },
        )
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            assert m["species_for_sbml"]["A"] == 10
        assert "species_for_sbml" in m

    def test_items(self):
        m = CompiledModel(species_for_sbml={"A": 1})
        items = dict(m.items())
        assert "species_for_sbml" in items
        assert items["species_for_sbml"] == {"A": 1}


class TestCompilerResult:
    def test_to_compiled_model(self):
        cr = CompilerResult(
            species_for_sbml={"A": 10},
            reactions_for_sbml={},
            parameters_for_sbml={},
        )
        cm = cr.to_compiled_model(
            species_not_mapped={"A": 10},
            mappings={"A": ["A"]},
        )
        assert isinstance(cm, CompiledModel)
        assert cm.species_for_sbml == {"A": 10}
        assert cm.species_not_mapped == {"A": 10}
        assert cm.mappings == {"A": ["A"]}

    def test_backward_compat_alias(self):
        cr = CompilerResult()
        cm = cr.to_compiled_model_dict(
            species_not_mapped={},
            mappings={},
        )
        assert isinstance(cm, CompiledModel)
