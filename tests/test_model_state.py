"""Regressions for model ownership, compilation, and execution boundaries."""

import math

import pytest

from mobspy import BaseSpecies, Simulation, Zero, u


def decay():
    A = BaseSpecies(["A"])
    A >> Zero @ (1 / u.minute)
    A(100)
    sim = Simulation(A)
    sim.duration = 1 * u.minute
    sim.volume = 0.2 * u.mL
    sim.plot_data = False
    sim.level = 0
    return A, sim


def test_recompilation_preserves_units_and_results():
    _, sim = decay()
    sim.compile(verbose=False)
    first = sim.generate_sbml()
    sim.run()
    expected = 500 / math.e
    assert sim.fres["A"][-1] == pytest.approx(expected, rel=1e-5)
    sim.compile(verbose=False)
    assert sim.generate_sbml() == first
    sim.run()
    assert sim.fres["A"][-1] == pytest.approx(expected, rel=1e-5)
    assert sim.parameters["duration"].units == u.minute.units
    assert sim.parameters["volume"].units == u.mL.units


@pytest.mark.parametrize("use_name", [False, True])
def test_update_concentration_updates_canonical_model_and_execution(use_name):
    A, sim = decay()
    sim.compile(verbose=False)
    sim.update_model(["A" if use_name else A, 500 / u.mL])
    assert sim.compiled_model.species["A"] == pytest.approx(100)
    sim.run(output_concentration=False)
    assert sim.fres["A"][0] == pytest.approx(100)
    assert sim.fres["A"][-1] == pytest.approx(100 / math.e, rel=1e-5)


def test_update_species_by_characteristic_name():
    A = BaseSpecies(["A"])
    A.alive >> A.dead @ 1
    sim = Simulation(A)
    sim.compile(verbose=False)
    sim.update_model(["A.alive", 50], ["A_dot_dead", 20])
    assert sim.compiled_model.species["A.alive"] == 50
    assert sim.compiled_model.species["A.dead"] == 20


def test_composition_repeated_exports_and_runs_do_not_modify_children():
    A = BaseSpecies(["A"])
    A >> Zero @ 1
    A(100)
    first, second = Simulation(A), Simulation(A)
    first.duration = second.duration = 1
    first.plot_data = second.plot_data = False
    composition = first + second
    assert len(composition.generate_sbml()) == 2
    assert len(composition.generate_sbml()) == 2
    for _ in range(2):
        composition.run()
        assert composition.fres["Time"][-1] == pytest.approx(2)
        assert composition.fres["A"][-1] == pytest.approx(100 / math.e**2, rel=1e-5)
    first.run()
    assert first.fres["Time"][-1] == pytest.approx(1)
    assert first.fres["A"][-1] == pytest.approx(100 / math.e, rel=1e-5)


def test_deleting_a_simulation_preserves_shared_declarations():
    A, sim = decay()
    sim.delete()
    another = Simulation(A)
    another.compile(verbose=False)
    assert len(another.compiled_model.reactions) == 1
    assert another.compiled_model.species["A"] == 100


def test_model_owns_only_its_declarations_and_survives_session_reset():
    from mobspy import Model
    from mobspy.dsl.session_context import reset_session

    A, B = BaseSpecies(["A", "B"])
    A >> Zero @ 1
    B >> Zero @ 2
    A(10)
    B(20)
    model = Model(A)
    reset_session()
    sim = Simulation(model)
    sim.compile(verbose=False)
    assert sim.compiled_model.species == {"A": 10}
    assert len(sim.compiled_model.reactions) == 1
    sim.delete()
    other = Simulation(model)
    other.compile(verbose=False)
    assert other.compiled_model.species == {"A": 10}


def test_session_does_not_retain_unreachable_models():
    import gc
    import weakref

    from mobspy.dsl.declarations import get_registry

    A, sim = decay()
    reference = weakref.ref(A)
    del A, sim
    gc.collect()
    assert reference() is None
    assert get_registry().reactions == []
    assert get_registry().counts == []


def test_composition_converts_stage_units_and_preserves_amounts():
    A = BaseSpecies(["A"])
    A >> Zero @ (1 / u.minute)
    A(100)
    first, second = Simulation(A), Simulation(A)
    first.duration, first.volume = 1 * u.minute, 1 * u.mL
    second.duration, second.volume = 60 * u.second, 0.002 * u.L
    combined = first + second
    combined.run(plot_data=False, unit_x=u.minute, unit_y=1 / u.mL)
    assert combined.fres["Time"][-1] == pytest.approx(2)
    assert combined.fres["A"][-1] == pytest.approx(50 / math.e**2, rel=1e-5)


def test_custom_backend_receives_compiled_models_and_resolved_settings():
    from mobspy.types import ExecutionPlan

    class Backend:
        def generate_model(self, model):
            return f"species={sorted(model.species)}"

        def run(self, plan, jobs=-1):
            assert isinstance(plan, ExecutionPlan)
            assert plan.settings[0].duration == 1
            assert plan.settings[0].volume == pytest.approx(0.2)
            assert plan.models[0][0].species["A"] == 100
            return [[{"Time": [0, 1], "A": [500, 250]}]]

    A, _ = decay()
    sim = Simulation(A, backend=Backend())
    sim.duration, sim.volume = 1 * u.minute, 0.2 * u.mL
    sim.run(plot_data=False)
    assert sim.fres["A"][-1] == 250
    assert sim.generate_sbml() == ["species=['A']"]


@pytest.mark.parametrize("context_manager", [False, True])
def test_event_times_keep_units_across_recompilation(context_manager):
    A, sim = decay()
    if context_manager:
        with sim.event_time(30 * u.second):
            A(0)
    else:
        sim.at(30 * u.second, {A: 0})
    for _ in range(2):
        sim.compile(verbose=False)
        sim.run()
        assert sim.fres["A"][-1] == pytest.approx(0, abs=1e-8)


def test_compiled_parameter_updates_do_not_change_a_shared_model():
    from mobspy import Model, ModelParameters

    A = BaseSpecies(["A"])
    k = ModelParameters(1)
    A >> Zero @ k
    model = Model(A)
    first = Simulation(model)
    first.compile(verbose=False)
    first.update_model([k, 2])
    second = Simulation(model)
    second.compile(verbose=False)
    assert first.compiled_model.parameters["k"][0] == 2
    assert second.compiled_model.parameters["k"][0] == 1


def test_exports_recompile_after_direct_configuration_changes():
    _, sim = decay()
    before = sim.generate_sbml()
    sim.parameters["volume"] = 2 * u.mL
    assert sim.generate_sbml() != before
    assert sim.compiled_model.unit_context.resolved_volume_magnitude == 2


def test_antimony_exports_every_stage_without_changing_sbml():
    _, first = decay()
    _, second = decay()
    composed = first + second
    before = composed.generate_sbml()
    assert len(composed.generate_antimony()) == 2
    assert composed.generate_sbml() == before


def test_numeric_duration_clears_previous_end_condition():
    A, sim = decay()
    sim.duration = A < 50
    sim.run()
    sim.duration = 2 * u.minute
    sim.run()
    assert sim.fres["Time"][-1] == pytest.approx(2)
    assert sim.fres["A"][-1] == pytest.approx(500 / math.e**2, rel=1e-5)


def test_deleted_simulation_releases_its_model():
    import gc
    import weakref

    A, sim = decay()
    reference = weakref.ref(A)
    sim.compile(verbose=False)
    sim.delete()
    del A
    gc.collect()
    assert reference() is None


def test_species_reset_removes_owned_reactions_after_session_reset():
    from mobspy.dsl.session_context import reset_session

    A, _ = decay()
    reset_session()
    A.reset_reactions()
    assert not A.get_reactions()
    sim = Simulation(A)
    sim.compile(verbose=False)
    assert not sim.compiled_model.reactions


@pytest.mark.parametrize("composition", [False, True])
def test_save_data_respects_filename_and_owns_results(tmp_path, composition):
    import json

    _, first = decay()
    sim = first + first if composition else first
    sim.output_dir = str(tmp_path / "new_directory")
    sim.output_file = "decay"
    result = sim.run(save_data=True)
    assert result is sim.results
    filename = tmp_path / "new_directory" / "decay.json"
    assert json.loads(filename.read_text()) == sim.results.to_dict()
    if composition:
        assert first.__dict__["results"] == {}
        assert first.plot_parameters.get("unit_x") is None


def test_composed_sbml_matches_sequential_run():
    import basico

    A, first = decay()
    first.duration = 1 * u.minute
    A.reset_reactions()
    A >> Zero @ (2 / u.minute)
    second = Simulation(A)
    second.duration = 2 * u.minute
    second.volume = 0.2 * u.mL
    chain = first + second
    chain.run(plot_data=False)
    expected = chain.fres["A"][-1]
    assert expected == pytest.approx(500 / math.e**5, rel=1e-5)
    sbml = chain.generate_sbml(compose=True)[0]
    model = basico.model_io.load_model_from_string(sbml)
    try:
        actual = basico.run_time_course(3, model=model, r_tol=1e-9, a_tol=1e-11)
        assert actual["A"].iloc[-1] == pytest.approx(expected, rel=1e-4)
    finally:
        basico.model_io.remove_datamodel(model)


def test_composed_export_rejects_mixed_units():
    from mobspy import SBMLError

    _, first = decay()
    _, second = decay()
    second.duration = 60 * u.second
    with pytest.raises(SBMLError, match="identical units and volumes"):
        (first + second).generate_sbml(compose=True)


def test_unit_parameter_updates_and_zero_values_use_model_time():
    from mobspy import ModelParameters

    A = BaseSpecies(["A"])
    k = ModelParameters(1 / u.minute)
    A >> Zero @ k
    A(100)
    sim = Simulation(A)
    sim.duration, sim.volume = 1 * u.minute, 0.2 * u.mL
    sim.compile(verbose=False)
    sim.update_model([k, [0 / u.minute, 2 / u.minute]])
    sim.run(plot_data=False)
    assert sim.results["A"][0][-1] == pytest.approx(500)
    assert sim.results["A"][1][-1] == pytest.approx(500 / math.e**2, rel=1e-5)
    assert sim.results.ts_model_parameters == [{"k": 0}, {"k": 2}]


def test_shared_sweeps_convert_independently_for_each_stage():
    from mobspy import ModelParameters

    A = BaseSpecies(["A"])
    k = ModelParameters([1 / u.minute, 2 / u.minute])
    A >> Zero @ k
    A(100)
    first, second = Simulation(A), Simulation(A)
    first.duration, second.duration = 1 * u.minute, 60 * u.second
    first.volume = second.volume = 0.2 * u.mL
    chain = first + second
    chain.run(plot_data=False)
    assert len(chain.results) == 2
    assert chain.results["A"][0][-1] == pytest.approx(500 / math.e**2, rel=1e-5)
    assert chain.results["A"][1][-1] == pytest.approx(500 / math.e**4, rel=1e-5)
    assert chain.results.ts_model_parameters == [{"k": 1}, {"k": 2}]


def test_concentration_parameter_sweeps_bind_initial_amounts():
    from mobspy import ModelParameters

    A = BaseSpecies(["A"])
    initial = ModelParameters([500 / u.mL, 1000 / u.mL])
    A >> Zero @ (1 / u.minute)
    A(initial)
    sim = Simulation(A)
    sim.duration, sim.volume = 1 * u.minute, 0.2 * u.mL
    sim.run(plot_data=False, output_concentration=False)
    assert sim.results["A"][0][0] == pytest.approx(100)
    assert sim.results["A"][1][0] == pytest.approx(200)
    assert sim.results["A"][1][-1] == pytest.approx(200 / math.e, rel=1e-5)


def test_second_order_unit_parameter_matches_analytic_decay():
    from mobspy import ModelParameters

    A = BaseSpecies(["A"])
    k = ModelParameters(1 * u.mL / u.minute)
    2 * A >> Zero @ k
    A(100)
    sim = Simulation(A)
    sim.duration, sim.volume = 0.01 * u.minute, 0.2 * u.mL
    sim.run(plot_data=False)
    assert sim.fres["A"][-1] == pytest.approx(500 / 11, rel=1e-5)


def test_replacing_an_end_condition_invalidates_compilation():
    A, sim = decay()
    sim.duration = A < 50
    sim.run()
    first_time = sim.fres["Time"][-1]
    sim.duration = A < 20
    sim.run()
    assert sim.fres["Time"][-1] > first_time
    assert sim.fres["A"][-1] <= 20
