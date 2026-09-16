# Simulation configuration

Assign parameters on a simulation or pass a configuration dictionary:

```python
S.duration = 10 * u.minute
S.volume = 0.2 * u.mL
S.run(plot_data=False, output_concentration=False)
```

`run()` keyword arguments update the user configuration. Compilation converts a
private copy into model units. The supplied quantities remain in `S.parameters`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `duration` | `60` | Duration in model time units, a time quantity, or a stop condition |
| `volume` | `1` | Compartment size; may carry volume units |
| `simulation_method` | `"deterministic"` | Solver method; `method` is an alias |
| `rate_type` | `None` | Infer deterministic or stochastic rate expansion |
| `repetitions` | `1` | Number of independent trajectories |
| `seeds` | `None` | One random seed per repetition |
| `jobs` | `-1` | Requested worker count |
| `step_size` | `None` | Optional output step size |
| `start_time` | `0` | Start time for output |
| `r_tol`, `a_tol` | `1e-8`, `1e-10` | Solver tolerances |
| `unit_x`, `unit_y` | `None` | Requested time and concentration/amount output units |
| `output_concentration` | `True` | Return concentrations; false requests amounts |
| `output_event` | `False` | Include event output points |
| `plot_data` | `True` | Plot after execution |
| `save_data` | `False` | Save results after execution |
| `output_dir` | `"outputs/"` | Directory for automatically saved results |
| `output_file` | `None` | Optional result filename |
| `level` | `2` | Logging level: 0 errors, 1 warnings, 2 information, 3 debug |

Use `S.plot_config` for plot settings and `S.set_from_json(path)` to load simulation
parameters. See `mobspy.simulation_config.SimulationConfig` for the complete
configuration definition.
