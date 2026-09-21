'''
Run the plotting examples headless -- a smoke test that a fresh install can
actually exercise the models end to end (build + data files + DataManager).

The examples double as visual verification: ``python examples/plot_<name>.py``.
'''
import os
import runpy
import pathlib

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

EX = pathlib.Path(__file__).resolve().parent.parent / "examples"

# examples known to run without network / private data (the local DuckDB store
# auto-fills what they need); the rest are work-in-progress plotting scripts.
RUNNABLE = ["plot_hltwim_fig12.py", "plot_hwm14_fig03.py", "plot_jvdm1_fig04.py",
            "plot_rocsat_fig03.py", "plot_sfqq_fig06.py", "plot_sfdd_fig07.py",
            "plot_sfpp_fig04.py", "plot_mm_fig02.py",
            "plot_eef_fig06.py", "plot_eej_fig02.py", "plot_survey_track.py",
            "plot_weimer05_fig02.py"]


@pytest.mark.parametrize("name", RUNNABLE)
def test_example_runs(name, monkeypatch):
    path = EX / name
    if not path.exists():
        pytest.skip("%s not present" % name)
    if "hltwim_fig12" in name:
        pytest.importorskip("apexpy")
    monkeypatch.setenv("MPLBACKEND", "Agg")
    monkeypatch.setattr(plt, "show", lambda *a, **k: None)
    monkeypatch.setattr("sys.argv", [str(path)])
    monkeypatch.chdir(EX.parent)
    try:
        runpy.run_path(str(path), run_name="__main__")
    finally:
        plt.close("all")
