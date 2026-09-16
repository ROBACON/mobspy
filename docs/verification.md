# Testing notebooks and examples

CI installs the built wheel and runs the full test suite on Python 3.11–3.14.
This includes all 16 tutorial notebooks and all 25 MobsPy example scripts in
`docs/example_models/`, including the application and journal subdirectories.
New notebooks and scripts in those directories are discovered automatically.

Each notebook runs every cell in a fresh Jupyter kernel using the same Python
environment as pytest, with Jupyter's inline plotting backend. Its working
directory is temporary, so exports cannot
depend on files left by another notebook. Errors fail the test, including cells
tagged as allowing errors. The inheritance, assignments, and plotting tutorials
also assert the results they demonstrate.

Scripts execute as written in separate Python processes, with copies of their
example assets in temporary directories. Matplotlib uses its headless `Agg`
backend, and Plotly uses its JSON renderer. Plot construction and file saving
remain enabled. The XOR example runs its full 50-by-50 parameter sweep; allow
several minutes for the complete suite.

## Run locally

From a checkout, install the test dependencies and package:

```bash
python -m pip install -e ".[test]"
python -m pytest tests/
```

Run only the notebooks and example scripts:

```bash
python -m pytest tests/test_notebooks.py tests/test_example_models.py
```

Jupyter must be allowed to open local sockets for communication with its kernels.
Executed notebooks, including output and failure tracebacks, are retained in
pytest's temporary directory. To explore the tutorials without pytest, install
`mobspy[examples]` in your notebook environment.

For release validation, CI builds a wheel, installs it with its `[test]` extra,
and runs `python -I -m pytest tests/ --import-mode=importlib --cov`. Isolated
Python processes and kernels prevent the source checkout from shadowing the
installed wheel.

## Other modeling tools

The `comparisons/` directory contains external reference models. A separate CI
job on Python 3.12 installs the versions in the `[comparisons]` extra and runs
BioCRNpyler, PySB, and BioNetGen. Use a separate virtual environment for these
tools because they have their own dependency constraints:

```bash
python -m pip install -e ".[comparisons]"
python docs/example_models/comparisons/comparison_bioCRNpyler.py
MPLBACKEND=Agg python docs/example_models/comparisons/hello_pysb.py
bionetgen run -i docs/example_models/comparisons/BioNetGen_multi_state.bngl -o /tmp/mobspy-bionetgen
```

The Kappa reference is an image, not an executable model. The MobsPy Kappa
comparison script is included in the regular example tests.

These checks establish that the examples execute, plots are constructed, exports
succeed, and tutorial assertions hold. They do not check interactive browser
controls or establish numerical equivalence between independent modeling tools.
