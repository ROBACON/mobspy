"""Execute every MobsPy example as shipped, including plots and full sweeps."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "docs" / "example_models"
MODELS = sorted(
    [
        *EXAMPLES_DIR.glob("*.py"),
        *(EXAMPLES_DIR / "application_models").glob("*.py"),
        *(EXAMPLES_DIR / "journal_models").glob("*.py"),
    ]
)


@pytest.mark.slow
@pytest.mark.parametrize(
    "model_path", MODELS, ids=[str(p.relative_to(EXAMPLES_DIR)) for p in MODELS]
)
def test_example_model(model_path: Path, tmp_path: Path) -> None:
    """Keep generated files isolated while exercising real plotting and saving."""
    examples = tmp_path / "examples"
    shutil.copytree(EXAMPLES_DIR, examples)
    script = examples / model_path.relative_to(EXAMPLES_DIR)
    env = {
        **os.environ,
        "MPLBACKEND": "Agg",
        "MPLCONFIGDIR": str(tmp_path / "matplotlib"),
        "PLOTLY_RENDERER": "json",
        "PYTHONHASHSEED": "0",
    }
    # The full 50-by-50 XOR parameter sweep takes about two minutes locally.
    timeout = 600 if model_path.name == "Logical_XOR_Gate.py" else 120
    result = subprocess.run(
        [sys.executable, "-I", str(script)],
        cwd=script.parent,
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )
    assert result.returncode == 0, (
        f"{model_path.relative_to(EXAMPLES_DIR)} failed:\n"
        f"STDOUT:\n{result.stdout[-4000:]}\nSTDERR:\n{result.stderr[-4000:]}"
    )
