"""Run example model scripts to verify they execute without errors.

Each example runs in a subprocess to avoid global state contamination
between MobsPy models.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES_DIR = Path(__file__).parent.parent / "docs" / "example_models"
APP_DIR = EXAMPLES_DIR / "application_models"
JOURNAL_DIR = EXAMPLES_DIR / "journal_models"

# Models that require external deps (plotly, pysb, seaborn) or are too slow
SKIP_MODELS = {
    "AND_gate.py",  # requires plotly
    "NOR_gate.py",  # requires plotly
    "hello_pysb.py",  # requires pysb
    "PySB_Comparison.py",  # requires pysb
    "Logical_XOR_Gate.py",  # requires seaborn + very slow parametric sweep
    "Phage_Transmission_System.py",  # requires local plot_config_donor.json
}

# Models that only compile (no .run()), used for comparison with other tools
COMPILE_ONLY = {
    "donor_receptor.py",
    "Antimony_Comparison.py",
    "BioCRNpyler_1.py",
    "BioCRNpyler_2.py",
    "BioNetGen_Comparison.py",
    "comparision_bioCRNpyler.py",
    "For_The_Trees.py",
    "Kappa_comp_v1.py",
}


def _collect_models(directory: Path) -> list[Path]:
    """Collect .py model files from a directory, excluding skipped ones."""
    if not directory.exists():
        return []
    return sorted(p for p in directory.glob("*.py") if p.name not in SKIP_MODELS)


def _run_model(model_path: Path, timeout: int = 120) -> None:
    """Run a model script in a subprocess with plotting/saving disabled."""
    # Create a wrapper script that patches plotting/saving then runs the model
    wrapper = f"""
import sys, os
import matplotlib
matplotlib.use('Agg')
os.environ['MPLBACKEND'] = 'Agg'

# Patch Simulation to disable plotting and saving by default
import mobspy.simulation as _sim_mod
_orig_init = _sim_mod.Simulation.__init__
def _patched_init(self, *args, **kwargs):
    _orig_init(self, *args, **kwargs)
    self.plot_data = False
    self.save_data = False
_sim_mod.Simulation.__init__ = _patched_init

# Suppress plt.show()
import matplotlib.pyplot as plt
plt.show = lambda *a, **k: None

import runpy
runpy.run_path({str(model_path)!r}, run_name='__main__')
"""
    result = subprocess.run(
        [sys.executable, "-c", wrapper],
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=model_path.parent,
        check=False,
    )
    assert result.returncode == 0, (
        f"Model {model_path.name} failed:\n"
        f"STDOUT:\n{result.stdout[-500:] if result.stdout else '(empty)'}\n"
        f"STDERR:\n{result.stderr[-500:] if result.stderr else '(empty)'}"
    )


# Application models
APPLICATION_MODELS = _collect_models(APP_DIR)


@pytest.mark.slow
@pytest.mark.parametrize(
    "model_path",
    APPLICATION_MODELS,
    ids=[p.stem for p in APPLICATION_MODELS],
)
def test_application_model(model_path: Path) -> None:
    """Run an application example model."""
    _run_model(model_path)


# Journal models
JOURNAL_MODELS = _collect_models(JOURNAL_DIR)


@pytest.mark.slow
@pytest.mark.parametrize(
    "model_path",
    JOURNAL_MODELS,
    ids=[p.stem for p in JOURNAL_MODELS],
)
def test_journal_model(model_path: Path) -> None:
    """Run a journal example model."""
    _run_model(model_path)
