"""Complete examples using MobsPy's new syntax.

Each example is self-contained and demonstrates one or more new features.
Run with: python docs/example_models/new_syntax_examples.py
"""

from __future__ import annotations


def example_1_at_rate_syntax() -> None:
    """Basic model using @ rate syntax."""
    from mobspy import BaseSpecies, Simulation

    A, B = BaseSpecies(["A", "B"])
    A >> B @ 0.1  # @ instead of []
    B >> A @ 0.05
    A(100)

    S = Simulation(A | B)
    S.duration = 50
    S.compile(verbose=False)
    print(f"Example 1: {len(S._reactions_for_sbml)} reactions compiled")


def example_2_reversible() -> None:
    """Reversible reaction using tuple rate."""
    from mobspy import BaseSpecies, Simulation

    A, B = BaseSpecies(["A", "B"])
    A >> B @ (1.0, 0.5)  # tuple = (forward, reverse)
    A(100)

    S = Simulation(A | B)
    S.duration = 20
    S.compile(verbose=False)
    print(f"Example 2: {len(S._reactions_for_sbml)} reactions (reversible)")


def example_3_explicit_events() -> None:
    """Events using S.at() and S.when()."""
    from mobspy import BaseSpecies, Simulation

    A, B = BaseSpecies(["A", "B"])
    A >> B @ 0.1
    A(100)

    S = Simulation(A | B)
    S.duration = 100

    # Time event: at t=20, set A to 50
    S.at(20, {A: 50})

    # Condition event: when A drops below 10, set B to 200
    S.when(A <= 10, {B: 200})

    S.compile(verbose=False)
    print(f"Example 3: {len(S._events_for_sbml)} events compiled")


def example_4_rate_builder() -> None:
    """Programmatic rate construction with the builder API."""
    from mobspy import BaseSpecies, Simulation
    from mobspy.modules.rate_builder import param_ref, species_ref

    A, B, C = BaseSpecies(["A", "B", "C"])

    # Build a rate expression as an AST (no lambda needed)
    rate = param_ref("k") * species_ref(A) * species_ref(B)
    A + B >> C @ rate

    A(50)
    B(30)

    S = Simulation(A | B | C)
    S.duration = 10
    S.compile(verbose=False)

    rxn = next(iter(S._reactions_for_sbml.values()))
    print(f"Example 4: kinetics = {rxn.kinetics}")


def example_5_hill_function() -> None:
    """Hill function rate using the builder."""
    from mobspy import BaseSpecies, Simulation, Zero
    from mobspy.modules.rate_builder import hill

    S_species, P = BaseSpecies(["S", "P"])

    # Enzymatic: Vmax * S^2 / (Km^2 + S^2)
    rate = hill("S", vmax=10, km=5, n=2)
    S_species >> P @ rate

    S_species(100)

    sim = Simulation(S_species | P)
    sim.duration = 20
    sim.compile(verbose=False)

    rxn = next(iter(sim._reactions_for_sbml.values()))
    print(f"Example 5: Hill kinetics = {rxn.kinetics}")


def example_6_conditional_rate() -> None:
    """Conditional rate using where()."""
    from mobspy import BaseSpecies, Simulation
    from mobspy.modules.rate_builder import literal, where

    A, B = BaseSpecies(["A", "B"])

    # Different rate depending on condition
    rate = where("A > 50", literal(0.1), literal(0.01))
    A >> B @ rate

    A(100)

    S = Simulation(A | B)
    S.duration = 50
    S.compile(verbose=False)

    rxn = next(iter(S._reactions_for_sbml.values()))
    print(f"Example 6: piecewise kinetics = {rxn.kinetics}")


def example_7_named_species() -> None:
    """Multi-axis species using .named()."""
    from mobspy import BaseSpecies, Simulation

    Color, Size = BaseSpecies(["Color", "Size"])
    Color.red
    Color.blue
    Size.big
    Size.small

    # Explicit naming (no sys._getframe)
    Widget = (Color * Size).named("Widget")
    Widget(10)
    Color.red >> Color.blue @ 0.1

    S = Simulation(Widget)
    S.compile(verbose=False)
    print(f"Example 7: species = {sorted(S._species_for_sbml.keys())}")


def example_8_standalone_compilation() -> None:
    """Compile and generate SBML without a Simulation instance."""
    from mobspy import BaseSpecies, compile_model
    from mobspy.execution import generate_sbml_from_compiled
    from mobspy.modules.list_species import List_Species
    from mobspy.modules.species_utils import create_orthogonal_vector_structure

    A, B = BaseSpecies(["A", "B"])
    A >> B @ 1.0
    A(10)

    model = List_Species([A, B])
    ovs = create_orthogonal_vector_structure(model)
    counts = [
        {
            "object": spe,
            "characteristics": c["characteristics"],
            "quantity": c["quantity"],
        }
        for spe in model
        for c in spe.get_quantities()
    ]
    reactions = set()
    for spe in model:
        for ref in spe.get_references():
            reactions.update(ref.get_reactions())

    result = compile_model(model, reactions, counts, ovs, verbose=False)

    sbml = generate_sbml_from_compiled(result)
    print(f"Example 8: SBML generated ({len(sbml)} chars)")


if __name__ == "__main__":
    examples = [
        example_1_at_rate_syntax,
        example_2_reversible,
        example_3_explicit_events,
        example_4_rate_builder,
        example_5_hill_function,
        example_6_conditional_rate,
        example_7_named_species,
        example_8_standalone_compilation,
    ]
    for ex in examples:
        ex()
    print("\nAll examples passed.")
