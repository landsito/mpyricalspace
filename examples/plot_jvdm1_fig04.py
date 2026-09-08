"""
Equatorial quiet-time vertical drift -- Figures 4 & 5 of Alken (2009)
-------------------------------------------------------------------

Fig. 4  "(left) Scherliess and Fejer model output ... (right) JULIA Vertical
        Drift Model output using an EUVAC index of 80."
Fig. 5  "Comparison of Scherliess and Fejer model output (blue), JVDM model
        output (green) ... as a function of season ... around (left) 1100 and
        (right) 1500 local times."

The raw JULIA data / its running mean (paper's middle panel and red curve) are
dropped -- only the two models are shown, in the same format. EUVAC 80 -> F10.7 80
for both models; drifts at Jicamarca, 150 km.

Both models are driven through ``Predictor.Empirical`` (the public interface).
The Scherliess & Fejer *quiet* drift is obtained from ``run_scherliessfejer`` by
feeding a flat 130 nT AE history (so the storm terms vanish) and reading the
``qvdrift`` field, one call per local time.

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_jvdm1_fig04.py
"""
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

from mpyricalspace import Predictor

JLON = -76.859012                                      # Jicamarca
FLUX = 80.0                                             # EUVAC index of 80
YEAR = 2015
LT = np.arange(8.0, 16.01, 0.25)
DOY = np.arange(1, 366)
QUIET_AE = np.array([[datetime(YEAR - 1, 12, 1), 130.0],
                     [datetime(YEAR + 1, 2, 1), 130.0]], dtype=object)
month = lambda x, _: (datetime(2000, 1, 1) + timedelta(days=int(x) - 1)).strftime("%b")[0]


def sf_quiet(lt, doy):
    """Scherliess & Fejer quiet-time drift [m/s] on the lt x doy grid (via Predictor)."""
    m = Predictor.Empirical()
    m.set_location(0.0, JLON, 0.0)
    out = np.empty((np.size(lt), np.size(doy)))
    for i, t in enumerate(lt):
        m.set_time([datetime(YEAR, 1, 1) + timedelta(days=int(d) - 1) for d in doy])
        out[i] = m.run_scherliessfejer(indices=QUIET_AE, F107=FLUX, SLT=float(t),
                                       hour_resolution=True).qvdrift.data
    return out


def jvdm_quiet(lt, doy):
    """Alken JVDM 150-km drift [m/s] on the lt x doy grid (via Predictor)."""
    m = Predictor.Empirical()
    m.set_time([datetime(YEAR, 1, 1)])
    ds = m.run_jvdm1(f107=FLUX, f107a=FLUX, slt=np.atleast_1d(lt), doy=doy)
    return np.atleast_2d(ds.qvdrift_150km.data).reshape(np.size(lt), np.size(doy))


SF = sf_quiet(LT, DOY)
JV = jvdm_quiet(LT, DOY)

# --------------------------------------------------------- Figure 4  (lt x season)
fig4, axes = plt.subplots(1, 2, figsize=(8, 4), sharey=True)
fig4.suptitle("Equatorial quiet-time vertical drift\n(Alken 2009, Fig. 4;  EUVAC 80)")
for ax, W, name in ((axes[0], SF, "Scherliess & Fejer"), (axes[1], JV, "Alken JVDM")):
    m = ax.pcolormesh(LT, DOY, W.T, vmin=-10, vmax=25, cmap="inferno", shading="auto")
    ax.set_title(name)
    ax.set_xlabel("local time (h)")
    ax.set_xlim(8, 16)
    ax.yaxis.set_major_locator(ticker.FixedLocator(range(15, 366, 30)))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(month))
axes[0].set_ylabel("season (month)")
fig4.colorbar(m, ax=axes, label="vertical drift (m/s)", shrink=0.85)
fig4.subplots_adjust(top=0.82)

# --------------------------------------------------------- Figure 5  (season, fixed LT)
fig5, axs = plt.subplots(1, 2, figsize=(9, 3.8))
fig5.suptitle("Equatorial quiet-time vertical drift vs season\n(Alken 2009, Fig. 5;  EUVAC 80)")
for ax, lt, ylim in ((axs[0], 11.0, (6, 26)), (axs[1], 15.0, (0, 14))):
    j = int(np.argmin(np.abs(LT - lt)))
    ax.plot(DOY, SF[j], color="tab:blue", label="Scherliess & Fejer")
    ax.plot(DOY, JV[j], color="tab:green", label="Alken JVDM")
    ax.set_title("LT = %04d" % (lt * 100))
    ax.set_xlabel("season (month)")
    ax.set_ylabel("vertical drift (m/s)")
    ax.set_xlim(1, 365)
    ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(ticker.FixedLocator(range(15, 366, 30)))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(month))
    ax.legend(fontsize=8)
fig5.tight_layout(rect=(0, 0, 1, 0.9))

plt.show()
