# MobsPy Migration Guide

This guide covers the new syntax and APIs introduced in the codebase refactor, with before/after examples for each change.

## Reaction Rates: `@` as alternative to `[]`

The `@` operator is an alternative syntax for attaching rates. Both `[]` and `@` are fully supported.

```python
# Original syntax (still works, no deprecation)
A >> B[1.5]
A + B >> C[lambda r1, r2: k * r1 * r2]

# Alternative syntax
A >> B @ 1.5
A + B >> C @ (lambda r1, r2: k * r1 * r2)
```

Multi-product reactions require parentheses (same as before):

```python
A >> (B + C) @ 0.5
(A + B) >> (C + D) @ 1.0
```

## Reversible Reactions: tuple replaces `Rev`

```python
# Before (deprecated, emits DeprecationWarning)
Rev[A >> B][1.0, 0.5]

# After (preferred)
A >> B @ (1.0, 0.5)   # tuple = (forward_rate, reverse_rate)
```

## Events: `S.at()` and `S.when()` replace context managers

```python
# Before (still works, not deprecated)
with S.event_time(10):
    A(50)
    B(0)

with S.event_condition(A <= 20):
    B(100)

# After (preferred - no modal operator behavior)
S.at(10, {A: 50, B: 0})
S.when(A <= 20, {B: 100})

# With delay
S.when(A <= 10, {B: 500}, delay=5)

# With characteristics
S.at(10, {A.alive: 50, A.dead: 0})
```

## Species Names: explicit replaces frame introspection

```python
# Before (deprecated, uses sys._getframe)
A, B = BaseSpecies(2)
C = New(A)
Thing = Color * Size

# After (preferred - no frame introspection)
A, B = BaseSpecies(['A', 'B'])
C = New(A, ['C'])
Thing = (Color * Size).named("Thing")
```

## Rate Builder: programmatic rate construction

The rate builder produces AST nodes directly, without lambda replay:

```python
from mobspy.modules.rate_builder import species_ref, param_ref, literal, where, hill

# Simple mass-action
rate = param_ref("k") * species_ref(A) * species_ref(B)
A + B >> C @ rate

# Conditional (equivalent to lambda r: 0.5 if ... else 1.0)
rate = where("A_dot_alive > 0", literal(0.5), literal(1.0))
A >> B @ rate

# Method syntax
rate = literal(0.5).where("A_dot_alive > 0", 1.0)

# Hill function
rate = hill("S", "Vmax", "Km", n=2)  # Vmax * S^n / (Km^n + S^n)

# Math functions
from mobspy.modules.rate_builder import rate_log, rate_exp, rate_sqrt
rate = rate_exp(-param_ref("k") * species_ref(A))
```

## Standalone Functions: use without Simulation

```python
from mobspy import compile_model, generate_sbml_from_compiled, run_sbml, plot_results
from mobspy.model_generation import generate_sbml_strings, generate_antimony_strings
from mobspy.modules.list_species import List_Species
from mobspy.modules.species_utils import create_orthogonal_vector_structure

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
# Before (deprecated, emits DeprecationWarning)
model["species_for_sbml"]
model["reactions_for_sbml"]["r1"]

# After (preferred)
model.species_for_sbml
model.reactions_for_sbml["r1"]
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
