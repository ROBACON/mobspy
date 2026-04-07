# MobsPy Migration Guide

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
# - generate_model(compiled, model_context=None) -> str
# - run(model_strings, parameters, jobs=-1) -> list
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
