"""
Manoj, Maus, Luehr & Alken (2008) prompt-penetration transfer function -- Fig. 2
------------------------------------------------------------------------------

"TF response to the synthetic inputs of IEF Ey. (left) The input signal and
(right) the zonal electric field output."

How the model normally works
    ``run_ppeefm1(nativelypackaged_code=True)`` (the local RTEEF) reads the interplanetary
    electric field IEF Ey from bundled ACE solar-wind data between the requested
    start / end times and pushes it, sample by sample, through a fixed 5-tap IIR
    filter (coefficients in ``TF.COF``:
        a = [1, -0.6023, -0.6600, 0.4451, -0.0463]
        b = [0.0052, 0.0151, 0.0014, -0.0152, -0.0061])
    -- equivalent to ``scipy.signal.lfilter(b, a, IEF_Ey)`` -- to obtain the
    prompt-penetration equatorial zonal electric field, with a ~17 min group
    delay. It also adds a climatological baseline (the CHAMP EEF model); F10.7
    only enters there, never the filter.

What this script does
    ``run_ppeefm1(ief=...)`` hands the model a *synthetic* IEF Ey instead: the
    array overrides the internal ACE buffer, so the identical TF.COF filter runs
    on made-up step / pulse / triangle inputs. Only the filter (prompt) output is
    meaningful in this mode -- the climatology term is skipped (``qef`` = 0).

    The input is 5-min cadence and must cover one hour of filter priming before
    the start through the end:  ``len(ief) == (dn - d0) / 300 s + 13``.

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_mm_fig02.py
"""
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt

from mpyricalspace import Predictor

STEP_H = 5.0 / 60.0                                     # 5-min cadence, in hours
T = np.arange(0.0, 6.0 + STEP_H / 2, STEP_H)            # 0..6 h output grid
PRIME = 12                                              # 1 h of filter-priming samples before t = 0

D0 = datetime(2004, 6, 1, 12)                           # throwaway window -- only its span matters
DN = D0 + timedelta(hours=6)


def box(a, b):
    return lambda t: ((t >= a) & (t < b)).astype(float)


def tri(a, peak, b):
    return lambda t: np.interp(t, [a, peak, b], [0.0, 1.0, 0.0], left=0.0, right=0.0)


CASES = [
    ("sustained step (on at 1 h)",  "tab:red",  box(1.0, 1e3)),
    ("box pulse, 1-3 h",            "tab:cyan", box(1.0, 3.0)),
    ("narrow triangle at 2 h",      "k",        tri(1.95, 2.0, 2.05)),
    ("wide triangle, 1.5-4.5 h",    "tab:blue", tri(2.5, 3.5, 4.5)),
]

_M = Predictor.Empirical()
_M.set_period(D0, DN, freq=timedelta(minutes=5))


def prompt_ef(shape):
    """prompt-penetration zonal E-field [mV/m] for a synthetic IEF Ey ``shape(t)``."""
    ief = np.concatenate((np.zeros(PRIME), shape(T)))   # PRIME priming zeros, then the signal
    return _M.run_ppeefm1(nativelypackaged_code=True, ief=ief).ppef.data


fig, (axl, axr) = plt.subplots(1, 2, figsize=(9, 3.8))
fig.suptitle("Manoj et al. (2008) prompt-penetration transfer function\n"
             "(reproducing their Fig. 2)")
for label, color, shape in CASES:
    axl.plot(T, shape(T), color=color, lw=2.0, label=label)
    axr.plot(T, prompt_ef(shape), color=color, lw=2.0)

axl.set_ylabel("IEF Ey (mV/m)"); axl.set_ylim(-1, 2)
axr.set_ylabel("equatorial zonal EF (mV/m)"); axr.set_ylim(-0.1, 0.1)
for ax in (axl, axr):
    ax.axhline(0, color="0.7", lw=0.6)
    ax.set_xlabel("time (hours)")
    ax.set_xlim(0, 6)
    ax.grid(True)
axl.legend(fontsize=7, loc="upper right")
fig.tight_layout(rect=(0, 0, 1, 0.88))

plt.show()
