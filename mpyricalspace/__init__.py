import importlib.metadata

from mpyricalspace._build import build_info, build_report

__version__ = importlib.metadata.version("mpyricalspace")

__all__ = ["__version__", "Predictor", "DataManager", "models", "survey",
           "build_info", "build_report"]

# Opt-in (MPYRICALSPACE_VSCODE_AUTOINSTALL=1), one-time, best-effort: install the
# companion VS Code extension when imported from a VS Code terminal (wheels have no
# post-install hook). Off by default; never blocks, never raises.
try:
    from mpyricalspace._vscode import maybe_autoinstall as _maybe_vscode
    _maybe_vscode(background=True)
except Exception:
    pass
