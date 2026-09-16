# MobsPy 3.0 Migration Guide

This guide covers the preferred syntax and APIs, with before/after examples for each change.

## Reaction Rates: `@` replaces `[]`

The `@` operator is the preferred syntax for attaching rates. `[]` still works but emits a `DeprecationWarning`.

```python
# Deprecated (emits DeprecationWarning)
A >> B[1.5]
A + B >> C[lambda r1, r2: k * r1 * r2]

# Preferred
A >> B @ 1.5
A + B >> C @ (lambda r1, r2: k * r1 * r2)
```

## Reversible Reactions: tuple replaces `Rev`

```python
# Deprecated (emits DeprecationWarning)
Rev[A >> B][1.0, 0.5]

# Preferred
A >> B @ (1.0, 0.5)   # tuple = (forward_rate, reverse_rate)
```

## Events: `S.at()` and `S.when()` replace context managers

```python
# Context manager style (still works, not deprecated)
with S.event_time(10):
    A(50)
    B(0)

with S.event_condition(A <= 20):
    B(100)

# Preferred -- explicit, no modal behavior
S.at(10, {A: 50, B: 0})
S.when(A <= 20, {B: 100})

# With delay
S.when(A <= 10, {B: 500}, delay=5)

# With characteristics
S.at(10, {A.alive: 50, A.dead: 0})
```

## Species Names: explicit replaces frame introspection

```python
# Uses sys._getframe (fragile in some environments)
A, B = BaseSpecies(2)
C = New(A)
Thing = Color * Size

# Preferred -- explicit names
A, B = BaseSpecies(['A', 'B'])
C = New(A, ['C'])
Thing = (Color * Size).named("Thing")
```

## Simulation Backend: pluggable via `backend` parameter

```python
from mobspy.sbml.backend import SBMLBackend

# Default (SBML/COPASI) -- same as omitting backend
S = Simulation(A | B, backend=SBMLBackend())

# Custom backends must implement the SimulationBackend protocol:
# - generate_model(model: ConcreteModel) -> str
# - run(plan: ExecutionPlan, jobs: int = -1) -> BackendResults
```

## Non-reactant species in rates

Any species can be referenced directly inside lambda rate functions, even if it is not a reactant:

```python
# Protein represses Gene (Protein is not a reactant in this reaction)
Gene >> Zero @ (lambda gene: 0.01 * gene * Protein ** 2 / (100 ** 2 + Protein ** 2))

# Or use the hill() convenience function:
Gene >> Zero @ hill(Protein, vmax=0.01, km=100, n=2, repression=True)

# Conditional rate
A >> B @ where("A > 50", 0.5, 1.0)
```

## Standalone Functions: use without Simulation

```python
from mobspy import compile_model, generate_sbml_from_compiled, run_sbml, plot_results
from mobspy.model_generation import generate_sbml_strings, generate_antimony_strings
from mobspy.dsl.list_species import List_Species
from mobspy.dsl.species_utils import create_orthogonal_vector_structure

# Build model
A, B = BaseSpecies(['A', 'B'])
A >> B @ 1.0
A(10)

# Compile without Simulation
model = List_Species([A, B])
ovs = create_orthogonal_vector_structure(model)
counts = [
    {"object": spe, "characteristics": c["characteristics"], "quantity": c["quantity"]}
    for spe in model for c in spe.get_quantities()
]
reactions = set()
for spe in model:
    for ref in spe.get_references():
        reactions.update(ref.get_reactions())

result = compile_model(model, reactions, counts, ovs, verbose=False)

# Generate SBML
sbml = generate_sbml_from_compiled(result)
```

## Dataclass Attribute Access: replaces dict-style

```python
# Deprecated (emits DeprecationWarning)
model["species_for_sbml"]
model["reactions_for_sbml"]["r1"]

# Preferred
model.species_for_sbml
model.reactions_for_sbml["r1"]
```

## Plot Configuration

```python
# PlotConfig is a dict subclass with attribute access
S.plot_config.xlabel = "Time (s)"
S.plot_config.ylabel = "Count"
S.plot_config["simulation_method"] = "stochastic"  # dict-style also works
```

## Structured Types

```python
from mobspy.types import ConcreteSpeciesId, CharacteristicQuery
from mobspy.constants import CharacteristicMarker

# ConcreteSpeciesId replaces _dot_ string surgery
sid = ConcreteSpeciesId("A", ("alive", "red"))
sid.to_sbml_id()    # "A_dot_alive_dot_red"
sid.to_display()    # "A.alive.red"
ConcreteSpeciesId.from_sbml_id("A_dot_alive")  # round-trip

# CharacteristicMarker enum replaces magic strings
CharacteristicMarker.ALL   # "all$" but type-safe
CharacteristicMarker.STD   # "std$"
CharacteristicMarker.NOT   # "not$"
```

## Semantic Testing

```python
from tests.conftest import CompiledModelAssertions

S.compile(verbose=False)
m = CompiledModelAssertions(S)
m.has_species("A", "B")
m.species_count("A", 10)
m.has_n_reactions(2)
m.has_reaction_involving("A")
m.kinetics_contains("1.5")
```


## Model ownership and cleanup

Capture a reusable model after declaring species, reactions, and initial counts:

```python
from mobspy import BaseSpecies, Model, Simulation, Zero

A = BaseSpecies(["A"])
A >> Zero @ 1
A(100)
model = Model(A)
first = Simulation(model)
second = Simulation(model)
first.delete()
second.run(plot_data=False)
```

`Model` captures reaction membership and initial counts. Define species
characteristics before capturing it. Simulations can compile a captured model in
another thread without consulting that thread's DSL registry. Configuration,
compiled parameter updates, and results belong to each simulation.

`delete()` releases only the receiving simulation's owned state. To remove a
reaction from future models, use the species' `reset_reactions()` method. A deleted
simulation cannot run again; create another simulation from the retained Model.
Deleting a composition clears its results and plan, leaving its children usable.

## Configuration and compiled output

Compilation keeps user-provided units in `S.parameters`. For example,
`S.duration = 1 * u.minute` remains a minute after compilation. Runtime conversion
happens on a separate configuration copy. Recompilation is safe, including after
changing volume or duration. Output time defaults to the supplied duration unit.

`S.compiled_model` exposes the authoritative `ConcreteModel`. Use
`S.update_model([A, count])` or `S.update_model(["A.alive", count])` for updates.
Parameter updates on one simulation do not change another simulation built from
the same Model. Legacy internal SBML attributes are derived views.

Compositions execute their stages without appending to the first child's internal
model list. Repeated exports and runs therefore preserve the same stage count.
Compositions also own their plot configuration; editing it does not change a child.
Stages may use different time and volume units; output uses the first stage's
unit system unless `unit_x` or `unit_y` requests another.

## Module paths

The supported public API is exported from `mobspy` and enumerated in
`mobspy.__all__`. Prefer those imports. Internal paths changed in 3.0:

| Previous location | Current location |
| --- | --- |
| `mobspy.modules.meta_class` | `mobspy.dsl.meta_class` |
| `mobspy.modules.compiler` | `mobspy.compiler.compiler` |
| `mobspy.modules.mobspy_expressions` | `mobspy.expressions.evaluation` |
| `mobspy.modules.unit_handler` | `mobspy.units.handler` |
| `mobspy.data_handler.time_series_object` | `mobspy.results.time_series` |
| `mobspy.parameter_scripts` | `mobspy.params` |
| `mobspy.plot_scripts` | `mobspy.plot` |
| `mobspy.sbml_simulator` | `mobspy.sbml` |

These module moves can include renamed symbols. Import documented public symbols
from the package root; consult the API documentation for internal helpers.
