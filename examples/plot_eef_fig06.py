"""
CHAMP equatorial zonal electric field climatology -- Figure 6 of
Alken & Maus (2010)
------------------------------------------------------------------

"[The CHAMP EEF model] electric field as a function of local time and longitude
for the four seasons."  (the paper's caption mislabels it "Fig. 1"; it is Fig. 6)

``run_eef`` is a function of longitude, local time, season (day of year), solar
flux and lunar local time.  Here it is evaluated on a longitude x local-time grid
for one representative day per season, at a fixed moderate flux, and **averaged
over lunar phase** (the lunar tide is a small ~0.1 mV/m modulation on top of the
solar-driven pattern shown here).  The model is only defined for 07-17 h local
time.

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_eef_fig06.py
"""
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from mpyricalspace import Predictor

FLUX = 110.0                                            # moderate solar flux
LON = np.arange(-180.0, 180.1, 5.0)
LT = np.arange(7.05, 17.001, 0.2)                       # inside the model's 7-17 h window
TAUS = np.linspace(0.0, 24.833, 8, endpoint=False)      # lunar local times to average over

SEASONS = [("March Equinox", 80), ("June Solstice", 172),
           ("September Equinox", 266), ("December Solstice", 355)]

_M = Predictor.Empirical()
_M.set_time([datetime(2005, 1, 1)])


def eef_season(doy):
    """lunar-phase-averaged zonal E-field [mV/m] on the (LON, LT) grid for one season."""
    acc = np.zeros((LON.size, LT.size))
    for tau in TAUS:
        acc += _M.run_eef(lon=LON, flux=FLUX, slts=LT, doys=[doy], lunars=float(tau)).eef.data
    return acc / len(TAUS)


fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True)
fig.suptitle("CHAMP equatorial zonal electric field\n(Alken & Maus 2010, Fig. 6;  F10.7 = %.0f)" % FLUX)
mesh = None
for ax, (name, doy) in zip(axes.flat, SEASONS):
    W = eef_season(doy)
    mesh = ax.pcolormesh(LON, LT, W.T, cmap="nipy_spectral", vmin=-0.8, vmax=1.0, shading="gouraud")
    ax.set_title(name)
    ax.set_xticks(range(-135, 136, 45))
    ax.set_yticks(range(8, 17, 2))
for ax in axes[:, 0]:
    ax.set_ylabel("local time (hours)")
for ax in axes[-1, :]:
    ax.set_xlabel("longitude (degrees)")

fig.subplots_adjust(right=0.88, top=0.86, hspace=0.22, wspace=0.08)
fig.colorbar(mesh, ax=axes, label="electric field (mV/m)", shrink=0.85,
             ticks=np.arange(-0.8, 1.01, 0.2))
plt.show()
