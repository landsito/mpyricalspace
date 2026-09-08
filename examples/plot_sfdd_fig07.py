"""
Scherliess & Fejer disturbance-dynamo drift -- Figs 4, 5, 7 & 8 of
Scherliess & Fejer (1997)
----------------------------------------------------------------

Fig. 4  "Local time variation of equatorial disturbance dynamo vertical drifts at
        the storm times shown in the top."  (an AE step, sampled at t1/t2/t3)
Fig. 5  "Comparison of the empirical disturbance dynamo drift pattern ... for an
        increase ... of about 400 nT over our quiet time level."  (tau = 9 h)
Fig. 7  "Local time variation of the long-term disturbance dynamo drifts for
        geomagnetically quiet short-term conditions."  (tau = 22-28 h)
Fig. 8  "Plots of the basic disturbance dynamo drift components ... for values of
        the AE index of 400 nT above our quiet time level."  (tau 1-6 / 1-12 / 22-28 h)

Driven through the public ``Predictor.Empirical`` interface: ``run_scherliessfejer``
is given a synthetic ``[datetime, AE]`` history (a step relative to the ~130 nT
quiet level) and its ``dynamo`` field is read.  The dynamo has no seasonal /
longitudinal dependence -- it is a pure functional of the AE time series.  The
1 h and 15 min resolutions share the same dynamo coefficients and AE windows, so
the resolution flag only matters for prompt penetration; the hourly model is used
here (``hour_resolution=True``).

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_sfdd_fig07.py
"""
from datetime import datetime, timedelta

from matplotlib import ticker
import numpy as np
import matplotlib.pyplot as plt

from mpyricalspace import Predictor

LT = np.arange(0.0, 24.001, 0.5)
QUIET, AE_D = 130.0, 400.0                              # quiet AE level, disturbance above it
ONSET = datetime(2015, 6, 1)                            # reference instant == "storm time 0"
_M = Predictor.Empirical()
_M.set_location(0.0, 0.0, 0.0)


def _grid(ax, major=5.0):
    """paper-style axes: inward ticks, major/minor y grid (comparison is direct)."""
    ax.tick_params(which="both", direction="in")
    ax.yaxis.set_major_locator(ticker.MultipleLocator(major))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(major / 2.0))
    ax.grid(True)


def _indices(intervals, peak):
    """[datetime, AE] step history: QUIET, then `peak` nT during each (a, b) storm-time-hour window."""
    idx = [[ONSET - timedelta(hours=200), QUIET]]
    for a, b in intervals:
        idx += [[ONSET + timedelta(hours=a), peak], [ONSET + timedelta(hours=b + 1), QUIET]]
    return np.array(sorted(idx, key=lambda r: r[0]), dtype=object)


def ae_of(intervals, storm_hours, peak=QUIET + AE_D):
    """the step's AE value (nT) at each storm-time hour -- for the AE panel."""
    h = np.asarray(storm_hours)
    out = np.full(h.shape, QUIET)
    for a, b in intervals:
        out[(h >= a) & (h <= b)] = peak
    return out


def dd(intervals, tau, peak=QUIET + AE_D):
    """disturbance-dynamo drift [m/s] over LT, `tau` h into a storm with this AE history."""
    idx = _indices(intervals, peak)
    _M.set_time([ONSET + timedelta(hours=tau)])
    return np.array([_M.run_scherliessfejer(indices=idx, hour_resolution=True,
                                            SLT=float(s)).dynamo.data[0] for s in LT])


# ------------------------------------------------------------------ Figure 4
STEP = [(0, 15)]                                        # 130 -> 400 nT at storm time 0, back at 15
TIMES = [("Storm Time t1", 6), ("Storm Time t2", 12), ("Storm Time t3", 21)]
fig4, ax = plt.subplots(4, 1, figsize=(5.5, 8))
fig4.suptitle("S&F disturbance dynamo -- an AE step\n(S&F 1997, Fig. 4)")
sthr = np.arange(-6, 25)
ax[0].step(sthr, ae_of(STEP, sthr, peak=400.0), where="post", color="k")
ax[0].set_ylabel("AE (nT)"); ax[0].set_xlabel("storm time (h)")
ax[0].set_xlim(-6, 24); ax[0].set_xticks(range(0, 25, 6)); _grid(ax[0], major=100.0)
for (label, tau) in TIMES:                              # t1 / t2 / t3 markers on the AE panel
    ax[0].axvline(tau, color="0.5", lw=0.8, ls="--")
    ax[0].annotate(label.split()[-1], xy=(tau, 400.0), xytext=(tau, 430.0),
                   ha="center", fontsize=8)
ax[0].set_ylim(90, 470)
for a, (label, tau) in zip(ax[1:], TIMES):
    a.axhline(0, color="0.7", lw=0.6)
    a.plot(LT, dd(STEP, tau, peak=400.0), "k-", lw=1.3)
    a.set_title(label, fontsize=9)
    a.set_xlim(0, 24); a.set_xticks(range(0, 25, 6))
    a.set_ylabel("drift (m/s)")
    _grid(a)
ax[-1].set_xlabel("local time (h)")
fig4.tight_layout(rect=(0, 0, 1, 0.92))

# ------------------------------------------------------------------ Figures 5 & 7
fig57, (a5, a7) = plt.subplots(1, 2, figsize=(9, 3.8))
fig57.suptitle("S&F disturbance dynamo drift patterns\n(S&F 1997, Figs. 5 & 7)")
a5.axhline(0, color="0.7", lw=0.6)
a5.plot(LT, dd([(0, 999)], 9), "k-", lw=1.3)
a5.set_title(r"$\tau$ = 9 h  (400 nT step)")
a7.axhline(0, color="0.7", lw=0.6)
a7.plot(LT, dd([(-32, -13)], 0), "k-", lw=1.3)          # active 22-28 h ago, quiet since
a7.set_title(r"$\tau$ = 22-28 h,  AE(1-12) $\approx$ 130 nT")
for a in (a5, a7):
    a.set_xlim(0, 24); a.set_xticks(range(0, 25, 6))
    a.set_xlabel("local time (h)"); a.set_ylabel("vertical drift (m/s)")
    _grid(a)
fig57.tight_layout(rect=(0, 0, 1, 0.84))

# ------------------------------------------------------------------ Figure 8
COMPONENTS = [(r"$\tau$ = 1-6 h", [(-6, 0)]),
              (r"$\tau$ = 1-12 h", [(-12, 0)]),
              (r"$\tau$ = 22-28 h,  AE(1-12) $\approx$ 130 nT", [(-32, -13)])]
fig8, ax8 = plt.subplots(3, 1, figsize=(5.5, 7), sharex=True)
fig8.suptitle("S&F basic disturbance-dynamo components,  AE$_d$ = 400 nT\n(S&F 1997, Fig. 8)")
for a, (label, intervals) in zip(ax8, COMPONENTS):
    a.axhline(0, color="0.7", lw=0.6)
    a.plot(LT, dd(intervals, 0), "k-", lw=1.3)
    a.set_title(label, fontsize=9)
    a.set_xlim(0, 24); a.set_xticks(range(0, 25, 6))
    a.set_ylabel("drift (m/s)")
    _grid(a)
ax8[-1].set_xlabel("local time (h)")
fig8.tight_layout(rect=(0, 0, 1, 0.9))

plt.show()
