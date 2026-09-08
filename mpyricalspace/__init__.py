import importlib.metadata

from mpyricalspace._build import build_info, build_report

__version__ = importlib.metadata.version("mpyricalspace")

__all__ = ["__version__", "Predictor", "DataManager", "models", "survey",
           "build_info", "build_report"]

# One-time, best-effort: offer the companion VS Code extension when imported from
# a VS Code terminal (wheels have no post-install hook). Never blocks, never
# raises; opt out with MPYRICALSPACE_NO_VSCODE=1.
try:
    from mpyricalspace._vscode import maybe_autoinstall as _maybe_vscode
    _maybe_vscode(background=True)
except Exception:
    pass
