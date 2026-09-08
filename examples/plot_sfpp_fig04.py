"""
Fejer & Scherliess prompt-penetration drift -- Figs 4 & 5 of
Fejer & Scherliess (1997)
---------------------------------------------------------

Fig. 4  "Empirical prompt penetration vertical drift patterns at three storm
        times following a step function increase in the AE index by 400 nT."
        (storm times t0 + 7.5 / 30 / 75 min)
Fig. 5  "Comparison of prompt penetration zonal electric fields obtained from our
        empirical model (solid curves) and from the Rice Convection Model ..."
        Only the empirical-model solid curves are reproduced here (initial time
        response and t + 60 min); the RCM curves are a separate model.

Driven through the public ``Predictor.Empirical`` interface: ``run_scherliessfejer``
is given a synthetic ``[datetime, AE]`` step (130 -> 130 + 400 nT at t0) and its
``prompt`` field is read.  The 15-min resolution is used (the default) -- its
three prompt terms are built from the AE differences over 0-15, 15-45 and
45-105 min, so the storm-time labels 7.5 / 30 / 75 min are the mid-points of
those windows and are selected by placing the evaluation time k * 15 min after t0.

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_sfpp_fig04.py
"""
from datetime import datetime, timedelta

from matplotlib import ticker
import numpy as np
import matplotlib.pyplot as plt

from mpyricalspace import Predictor

STEP_MIN = 15.0                                         # 15-min resolution
QUIET, DELTA_AE = 130.0, 400.0
LT = np.arange(0.5, 24.0, 0.5)
ONSET = datetime(2015, 6, 1)                            # t0 -- the AE step
AE_STEP = np.array([[ONSET - timedelta(hours=40), QUIET],
                    [ONSET, QUIET + DELTA_AE]], dtype=object)
_M = Predictor.Empirical()
_M.set_location(0.0, 0.0, 0.0)


def ae_value(minutes):
    """the step's AE (nT) at `minutes` of storm time -- for the AE panel."""
    return np.where(np.asarray(minutes) >= 0.0, QUIET + DELTA_AE, QUIET)


def pp(k):
    """prompt-penetration drift [m/s] over LT, k * 15 min after t0 (15-min model)."""
    _M.set_time([ONSET + timedelta(minutes=STEP_MIN * k)])
    return np.array([_M.run_scherliessfejer(indices=AE_STEP, SLT=float(s)).prompt.data[0]
                     for s in LT])


def _grid(ax, major=5.0):
    """paper-style axes: inward ticks, major/minor y grid so the comparison is direct."""
    ax.tick_params(which="both", direction="in")
    ax.yaxis.set_major_locator(ticker.MultipleLocator(major))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(major / 2.0))
    ax.grid(True)


def ae_panel(ax, marks):
    """the AE step, storm time in minutes, with dashed markers at `marks` (label, minutes)."""
    s = np.arange(-2, 9)
    ax.step(s * STEP_MIN, ae_value(s * STEP_MIN), where="post", color="k")
    for label, mins in marks:
        ax.axvline(mins, color="0.5", lw=0.8, ls="--")
        ax.annotate(label, xy=(mins, 530.0), xytext=(mins, 560.0), ha="center", fontsize=8)
    ax.set_ylabel("AE (nT)")
    ax.set_xlabel("storm time (min)")
    ax.set_xlim(-30, 120)
    ax.set_ylim(90, 620)
    _grid(ax, major=200.0)


# ------------------------------------------------------------------ Figure 4
#   dAEt_7P5 : step 0-15 min ago   -> k = 0   ("t0 + 7.5 min")
#   dAEt_30  : step 15-45 min ago  -> k = 2   ("t0 + 30 min")
#   dAEt_75  : step 45-105 min ago -> k = 5   ("t0 + 75 min")
TIMES4 = [(r"t$_0$ + 7.5 min", 0), (r"t$_0$ + 30 min", 2), (r"t$_0$ + 75 min", 5)]
fig4, ax4 = plt.subplots(4, 1, figsize=(5.5, 8))
fig4.suptitle(r"F&S prompt penetration -- an AE step,  $\Delta$AE = 400 nT"
              "\n(Fejer & Scherliess 1997, Fig. 4)")
ae_panel(ax4[0], [(lbl, k * STEP_MIN + 7.5) for lbl, k in TIMES4])
for a, (label, k) in zip(ax4[1:], TIMES4):
    a.axhline(0, color="0.7", lw=0.6)
    a.plot(LT, pp(k), "k-", lw=1.3)
    a.set_title(label, fontsize=9)
    a.set_xlim(0, 24)
    a.set_xticks(range(0, 25, 4))
    a.set_ylabel("vertical drift (m/s)")
    _grid(a)
ax4[-1].set_xlabel("local time (h)")
fig4.tight_layout(rect=(0, 0, 1, 0.93))

# ------------------------------------------------------------------ Figure 5
TIMES5 = [("initial time response", 0, 12.0, 8), ("t + 60 min", 4, 16.0, -14)]
fig5, (a_ae, a_v) = plt.subplots(2, 1, figsize=(5.5, 7),
                                gridspec_kw={"height_ratios": [1, 3]})
fig5.suptitle("F&S prompt-penetration drift, empirical model only"
              "\n(Fejer & Scherliess 1997, Fig. 5)")
ae_panel(a_ae, [(r"t$_\epsilon$", 0.0), (r"t$_{60}$", 60.0)])
a_v.axhline(0, color="0.7", lw=0.6)
for label, k, lt_lbl, dy in TIMES5:
    curve = pp(k)
    a_v.plot(LT, curve, "k-", lw=1.3)
    yv = curve[np.argmin(np.abs(LT - lt_lbl))]
    a_v.annotate(label, xy=(lt_lbl, yv), xytext=(0, dy), textcoords="offset points",
                 ha="center", fontsize=9)
a_v.set_xlim(0, 24)
a_v.set_xticks(range(0, 25, 4))
a_v.margins(y=0.16)
a_v.set_xlabel("local time (h)")
a_v.set_ylabel("vertical drift (m/s)")
_grid(a_v)
fig5.tight_layout(rect=(0, 0, 1, 0.92))

plt.show()
