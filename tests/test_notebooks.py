"""Execute the tutorials from a clean kernel against the installed package."""

from __future__ import annotations

import sys
from pathlib import Path

import nbformat
import pytest
from jupyter_client import KernelManager
from nbclient import NotebookClient

NOTEBOOK_DIR = (
    Path(__file__).resolve().parent.parent
    / "docs"
    / "example_models"
    / "tutorial_notebooks"
)
NOTEBOOKS = sorted(NOTEBOOK_DIR.glob("*.ipynb"))


@pytest.mark.slow
@pytest.mark.parametrize("path", NOTEBOOKS, ids=[p.stem for p in NOTEBOOKS])
def test_tutorial_notebook(path: Path, tmp_path: Path, monkeypatch) -> None:
    """Run every cell, including assertions, plotting, and model/file exports."""
    monkeypatch.setenv("MPLBACKEND", "module://matplotlib_inline.backend_inline")
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "matplotlib"))
    monkeypatch.setenv("IPYTHONDIR", str(tmp_path / "ipython"))
    monkeypatch.setenv("JUPYTER_RUNTIME_DIR", str(tmp_path / "jupyter"))
    monkeypatch.setenv("PLOTLY_RENDERER", "json")
    notebook = nbformat.read(path, as_version=4)
    # A globally registered python3 kernel may belong to a different environment.
    manager = KernelManager(kernel_name="python3")
    manager.kernel_spec.argv = [
        sys.executable,
        "-I",
        "-m",
        "ipykernel_launcher",
        "-f",
        "{connection_file}",
    ]
    client = NotebookClient(
        notebook,
        km=manager,
        timeout=120,
        force_raise_errors=True,
        resources={"metadata": {"path": str(tmp_path)}},
    )
    try:
        client.execute()
    finally:
        if manager.has_kernel:
            manager.shutdown_kernel(now=True)
        # Preserve executed cells and tracebacks in pytest's temporary directory.
        nbformat.write(notebook, tmp_path / path.name)
