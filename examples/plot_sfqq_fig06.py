"""
Scherliess & Fejer quiet-time equatorial vertical drift -- Figs 6, 7 & 8 of
Scherliess & Fejer (1999)
--------------------------------------------------------------------------

Fig. 6  "Comparison of the model predictions with average drift patterns ...
        for equinox, June solstice, and December solstice conditions for low,
        medium, and high solar flux periods."  (model curves only, at Jicamarca)
Fig. 7  "Solar cycle dependence of the empirical model drifts in the
        African-Indian (0-150 E), Pacific (150-210 E), Western American
        (210-300 E), and Brazilian (300-360 E) equatorial regions."
Fig. 8  "Empirical model results in six longitudinal sectors for low (dashed)
        and high (solid line) solar flux conditions."

This is the *quiet-time* half of the model.  It is driven through the public
``Predictor.Empirical`` interface: ``run_scherliessfejer`` is fed a flat 130 nT
AE history (so the storm terms vanish) and its ``qvdrift`` field is read.  The
model is a four-month climatology, so any day inside a season window gives the
same curve -- a representative day is used.  Longitudinal regions are sampled in
10 deg steps.

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_sfqq_fig06.py
"""
from datetime import datetime, timedelta

from matplotlib import ticker
import numpy as np
import matplotlib.pyplot as plt

from mpyricalspace import Predictor

JLON = -76.859012                                       # Jicamarca
YEAR = 2015
LT = np.arange(0.0, 24.001, 0.25)
SEASONS = [("MAR-APR / SEP-OCT", 90), ("MAY-AUG", 182), ("NOV-FEB", 15)]
QUIET_AE = np.array([[datetime(YEAR - 1, 12, 1), 130.0],
                     [datetime(YEAR, 6, 1), 130.0]], dtype=object)
_M = Predictor.Empirical()


def sf(lon, doy, sa):
    """quiet-time vertical drift [m/s] over ``LT`` for one longitude / season / solar flux.

    One ``run_scherliessfejer`` call: the evaluation times are picked so each maps
    exactly onto an ``LT`` value (t + lon/15 = LT), the AE history is flat 130 nT
    (storm terms -> 0) and F10.7 is fixed at ``sa``.
    """
    base = datetime(YEAR, 1, 1) + timedelta(days=int(doy) - 1)
    _M.set_time([base + timedelta(hours=float((t - lon / 15.0) % 24.0)) for t in LT])
    _M.set_location(0.0, lon, 0.0)
    return _M.run_scherliessfejer(indices=QUIET_AE, F107=float(sa), hour_resolution=True).qvdrift.data


# ---------------------------------------------------- Figure 6  (Jicamarca, 3 Sa x 3 seasons)
fig6, axes = plt.subplots(3, 3, figsize=(8, 6), sharex=True, sharey=True)
fig6.suptitle("Scherliess & Fejer quiet-time drift at Jicamarca\n(S&F 1999, Fig. 6)")
for r, sa in enumerate((200, 140, 80)):
    for c, (label, doy) in enumerate(SEASONS):
        ax = axes[r, c]
        ax.axhline(0, color="0.7", lw=0.6)
        ax.plot(LT, sf(JLON, doy, sa), "k-", lw=1.3)
        ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 6))
        if r == 0:
            ax.set_title(label, fontsize=10)
        if c == 2:
            ax.text(1.02, 0.5, "Sa = %d" % sa, transform=ax.transAxes, va="center", rotation=90)
axes[2, 1].set_xlabel("local time (h)")
axes[1, 0].set_ylabel("vertical drift (m/s)")
_=[_ax.tick_params(which='both',direction='in') for _ax in axes.flatten()]
_=[_ax.yaxis.set_major_locator(ticker.MultipleLocator(20)) for _ax in axes.flatten()]
_=[_ax.yaxis.set_minor_locator(ticker.MultipleLocator(10)) for _ax in axes.flatten()]
_=[_ax.grid(True) for _ax in axes.flatten()]
fig6.tight_layout(rect=(0, 0, 1, 0.92))

# ---------------------------------------------------- Figure 7  (solar-cycle dependence, 4 regions)
# Procedure (S&F 1999, Fig. 7): for each longitude sector the drift residual is
# derived by subtracting the model response for Sa = 100 from the one for
# Sa = 200 at 10 deg longitude steps, then averaging those residuals over the
# sector's longitude bin.
REGIONS = [("African-Indian (0-150 E)", (0, 150), "-"), ("Pacific (150-210 E)", (150, 210), ":"),
           ("Western American (210-300 E)", (210, 300), "--"), ("Brazilian (300-360 E)", (300, 360), (0, (6, 3)))]
SA_LO, SA_HI = 100.0, 200.0
fig7, axs = plt.subplots(3, 1, figsize=(6, 8), sharex=True)
fig7.suptitle("Scherliess & Fejer solar-cycle dependence\n(S&F 1999, Fig. 7)")
F7_SEASONS = [SEASONS[1], SEASONS[2], SEASONS[0]]        # June solstice, December solstice, equinox
for ax, (label, doy) in zip(axs, F7_SEASONS):
    ax.axhline(0, color="0.7", lw=0.6)
    for name, (lo, hi), ls in REGIONS:
        lons = np.arange(lo + 5.0, hi, 10.0)            # 10 deg steps across the sector
        resid = np.mean([sf(lon, doy, SA_HI) - sf(lon, doy, SA_LO) for lon in lons], axis=0)
        ax.plot(LT, resid, color="k", ls=ls, lw=1.1, label=name)
    ax.set_title(label, fontsize=10)
    ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 6))
    ax.set_ylabel(r"$V(S_a{=}200) - V(S_a{=}100)$  [m/s]")
axs[0].legend(fontsize=7, loc="upper left")
axs[-1].set_xlabel("local time (h)")
_=[_ax.tick_params(which='both',direction='in') for _ax in axs]
_=[_ax.yaxis.set_major_locator(ticker.MultipleLocator(20)) for _ax in axs]
_=[_ax.yaxis.set_minor_locator(ticker.MultipleLocator(10)) for _ax in axs]
_=[_ax.grid(True) for _ax in axs]
fig7.tight_layout(rect=(0, 0, 1, 0.93))

# ---------------------------------------------------- Figure 8  (6 sectors x 3 seasons, low/high Sa)
PHI = (0, 60, 120, 180, 240, 300)
fig8, ax8 = plt.subplots(len(PHI), len(SEASONS), figsize=(6.8, 8.8), sharex=True, sharey=True)
fig8.suptitle("Scherliess & Fejer quiet-time drift, low (dashed) / high (solid) Sa\n(S&F 1999, Fig. 8)")
for r, phi in enumerate(PHI):
    for c, (label, doy) in enumerate(SEASONS):
        ax = ax8[r, c]
        ax.axhline(0, color="0.7", lw=0.6)
        ax.plot(LT, sf(phi, doy, 90), "k--", lw=1.0)
        ax.plot(LT, sf(phi, doy, 180), "k-", lw=1.1)
        ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 8))
        ax.tick_params(labelsize=7)
        if r == 0:
            ax.set_title(label, fontsize=9)
        if c == 2:
            ax.text(1.03, 0.5, r"$\phi$ = %03d E" % phi, transform=ax.transAxes, va="center", rotation=90, fontsize=8)
ax8[-1, 1].set_xlabel("local time (h)")
ax8[len(PHI) // 2, 0].set_ylabel("vertical drift (m/s)")
_=[_ax.tick_params(which='both', direction='in') for _ax in ax8.flatten()]
_=[_ax.yaxis.set_major_locator(ticker.MultipleLocator(20)) for _ax in ax8.flatten()]
_=[_ax.yaxis.set_minor_locator(ticker.MultipleLocator(10)) for _ax in ax8.flatten()]
_=[_ax.grid(True) for _ax in ax8.flatten()]
fig8.tight_layout(rect=(0, 0, 1, 0.93))

plt.show()
