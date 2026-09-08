"""
ROCSAT-1 quiet-time equatorial vertical drift -- Figures 3 & 4 of Fejer et al. (2008)
-----------------------------------------------------------------------------------

Fig. 3  "Local time and longitude dependence of quiet time equatorial vertical
        drifts (positive upward) in eight longitudinal sectors for moderate
        solar flux conditions."
Fig. 4  "Local time, seasonal and longitudinal dependent equatorial quiet time
        vertical drift velocities for moderate solar flux conditions."

Moderate solar flux -> Sa = 150.  The seasons follow the paper: Nov-Feb, equinox
(months around both equinoxes), May-Aug -- the model is averaged over the days in
each window.

The model is driven through the public ``Predictor.Empirical`` interface:
``run_rocsat`` takes explicit ``lons``, ``slts`` and ``doys`` arrays and returns
an (lon, doy, slt) cube, which is averaged over ``doy``.

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_rocsat_fig03.py
"""
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from mpyricalspace import Predictor

FLUX = 150.0                                            # moderate solar flux (Sa)
LT = np.arange(0.0, 24.001, 0.5)
LON = np.arange(-180.0, 180.1, 5.0)
SECTORS = [-165, -120, -75, -30, 15, 60, 105, 150]      # Fig. 3 longitude panels

SEASONS = {
    "NOV-FEB": np.r_[np.arange(305, 366, 6), np.arange(1, 60, 6)],
    "EQUINOX": np.r_[np.arange(36, 97, 6), np.arange(219, 281, 6)],
    "MAY-AUG": np.arange(121, 244, 6),
}
_M = Predictor.Empirical()
_M.set_time([datetime(2015, 1, 1)])


def drift(lons, slts, doys):
    """ROCSAT-1 quiet drift [m/s] on the (lon, slt) grid, averaged over `doys`."""
    ds = _M.run_rocsat(f107s=FLUX, lons=np.asarray(lons, float),
                       doys=np.asarray(doys, int), slts=np.asarray(slts, float))
    return ds.qvdrift_rocsat.mean("doy").transpose("lon", "slt").data


D = {name: drift(LON, LT, doys) for name, doys in SEASONS.items()}

# ------------------------------------------------------- Figure 3  (8 sectors x 3 seasons)
fig3, axes = plt.subplots(len(SECTORS), len(SEASONS), figsize=(7, 8.5), sharex=True)
fig3.suptitle("ROCSAT-1 quiet-time vertical drift\n(Fejer et al. 2008, Fig. 3;  Sa=150, Kp<=3)")
for j, season in enumerate(SEASONS):
    for i, lon in enumerate(SECTORS):
        ax = axes[i, j]
        k = int(np.argmin(np.abs(LON - lon)))
        ax.axhline(0, color="0.7", lw=0.6)
        ax.plot(LT, D[season][k], "k.-", ms=2.5)
        ax.set_xlim(0, 24)
        ax.set_ylim(-30, 45)
        ax.set_yticks([-20, 0, 20, 40])
        ax.tick_params(labelsize=7)
        ax.text(0.97, 0.9, "%d°" % lon, transform=ax.transAxes, ha="right", va="top", fontsize=8)
        if i == 0:
            ax.set_title(season, fontsize=10)
        if i == len(SECTORS) - 1:
            ax.set_xlabel("local time (h)")
axes[len(SECTORS) // 2, 0].set_ylabel("vertical drift (m/s)")
fig3.tight_layout(rect=(0, 0, 1, 0.93))

# ------------------------------------------------------- Figure 4  (contours, 3 seasons)
LEVELS = np.arange(-45.0, 45.1, 10.0)
fig4, axs = plt.subplots(len(SEASONS), 1, figsize=(5.2, 8.5), sharex=True)
fig4.suptitle("ROCSAT-1 quiet-time vertical drift\n(Fejer et al. 2008, Fig. 4;  Sa=150, Kp<=3)")
cf = None
for ax, (season, W) in zip(axs, D.items()):
    cf = ax.contourf(LON, LT, W.T, levels=LEVELS, cmap="jet", extend="both")
    ax.contour(LON, LT, W.T, levels=LEVELS, colors="k", linewidths=0.3)
    ax.set_title(season)
    ax.set_ylabel("local time (h)")
    ax.set_ylim(0, 24)
    ax.set_yticks(range(0, 25, 6))
axs[-1].set_xlabel("longitude (°)")
axs[-1].set_xticks(range(-180, 181, 90))
fig4.subplots_adjust(right=0.85, top=0.88)
fig4.colorbar(cf, ax=axs, label="m/s", shrink=0.7)

plt.show()
