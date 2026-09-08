"""
HWM14 quiet-time thermospheric winds -- Figure 3 of Drob et al. (2015)
--------------------------------------------------------------------

"HWM14 zonal (top) and meridional (bottom) quiet time winds at 250 km altitude,
as a function of solar local time and geographic latitude. The models were
evaluated under December solstice (day of year 0), June solstice (day 180), and
combined equinox (average of days 90 and 270) conditions."

The wind is climatological -- the DWM disturbance part is subtracted (u - du) --
and averaged over geographic longitude (the paper shows only the HWM14 column;
its WINDII / HWM07 columns are separate models).

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_hwm14_fig03.py
"""
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt

from mpyricalspace import Predictor

YEAR = 2020
ALT = 250.0
AP = 4.0                                                # any value -- the disturbance is subtracted

LAT = np.arange(-60.0, 60.1, 4.0)
SLT = np.arange(0.0, 24.0, 1.0)
GLON = np.arange(0.0, 360.0, 30.0)                      # averaged over
LA, ST, LO = np.meshgrid(LAT, SLT, GLON, indexing="ij")
UT = np.mod(ST - LO / 15.0, 24.0)
N = LA.size

COLS = [("Dec solstice\n(day 0)", (0,)),
        ("Jun solstice\n(day 180)", (180,)),
        ("equinox\n(days 90 & 270)", (90, 270))]
LEVELS = np.arange(-160.0, 160.1, 20.0)


def quiet_wind(doys):
    """(eastward, northward) quiet wind [m/s] on the (LAT, SLT) grid, averaged over doys and GLON."""
    us, vs = [], []
    for d in doys:
        times = [datetime(YEAR, 1, 1) + timedelta(days=int(d)) + timedelta(hours=float(u))
                 for u in UT.ravel()]
        m = Predictor.Empirical()
        m.set_time(times)
        ds = m.run_hwm(ap=AP, lat=LA.ravel(), lon=LO.ravel(), alt=np.full(N, ALT), version=2014)
        us.append((ds.u.data - ds.du.data).reshape(LA.shape).mean(axis=2))
        vs.append((ds.v.data - ds.dv.data).reshape(LA.shape).mean(axis=2))
    return np.mean(us, axis=0), np.mean(vs, axis=0)


fig, axes = plt.subplots(2, 3, figsize=(10, 5.8), sharex=True, sharey=True)
fig.suptitle("HWM14 quiet-time winds at %.0f km\n(Drob et al. 2015, Fig. 3)" % ALT)
cf = None
for j, (label, doys) in enumerate(COLS):
    U, V = quiet_wind(doys)
    for i, (W, name) in enumerate(((U, "Eastward"), (V, "Northward"))):
        ax = axes[i, j]
        cf = ax.contourf(SLT, LAT, W, levels=LEVELS, cmap="RdBu_r", extend="both")
        cl = ax.contour(SLT, LAT, W, levels=LEVELS, colors="k", linewidths=0.4)
        ax.clabel(cl, cl.levels[::2], fmt="%d", fontsize=6)
        if i == 0:
            ax.set_title(label)
        if j == 0:
            ax.set_ylabel("%s\nGeo. latitude" % name)
        if i == 1:
            ax.set_xlabel("solar local time (h)")
        ax.set_xticks(range(0, 25, 6))
        ax.set_ylim(-60, 60)

fig.subplots_adjust(right=0.9, top=0.82)
fig.colorbar(cf, ax=axes, label="m/s", shrink=0.8)
plt.show()
