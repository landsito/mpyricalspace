"""
Equatorial electrojet climatology from CHAMP / Orsted / SAC-C -- Figs 2 & 5 of
Alken & Maus (2007)
--------------------------------------------------------------------------

Fig. 2  "Longitudinal dependence of EEJ for different seasons."
Fig. 5  "Seasonal dependence of EEJ for different longitudes."
Both: local time fixed at 10:30, EUVAC = 150, lunar local time 12.4 h (the
defaults of Alken's ``eej_plot`` tool that produced the figures).

``run_eej(version=2, model=...)`` selects which satellite the coefficients were
fit from -- ``"champ"`` (default), ``"oersted"`` or ``"sac-c"``.  SAC-C sampled a
single local time, so its J has no local-time dependence.  Amplitudes may differ
a little from the published figures (coefficient revisions / EUVAC convention);
the longitude and season structure is what is reproduced here.

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_eej_fig02.py
"""
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from mpyricalspace import Predictor

LT, EUVAC, LUNAR = 10.5, 150.0, 12.4                    # fixed, as in Alken & Maus (2007) Figs 2 & 5
SATS = [("CHAMP", "champ", "tab:red"),
        ("Ørsted", "oersted", "tab:blue"),
        ("SAC-C", "sac-c", "tab:green")]

_M = Predictor.Empirical()
_M.set_time([datetime(2005, 1, 1)])


def eej(model, lons, doys):
    """EEJ current density [A/m] for a satellite model over the (lons x doys) grid."""
    ds = _M.run_eej(version=2, model=model, lon=np.atleast_1d(lons), flux=EUVAC,
                    slts=[LT], doys=np.atleast_1d(doys), lunars=LUNAR)
    return np.atleast_1d(ds.eej.data).reshape(np.size(lons), np.size(doys))


# ----------------------------------------------------- Figure 2  (J vs longitude, 4 seasons)
LON = np.arange(-180.0, 180.1, 3.0)
SEASONS = [("March Equinox", 80), ("June Solstice", 172),
           ("September Equinox", 266), ("December Solstice", 355)]

fig2, ax2 = plt.subplots(4, 1, figsize=(5.5, 9), sharex=True)
fig2.suptitle("EEJ longitudinal dependence by season\n(Alken & Maus 2007, Fig. 2)")
for ax, (name, doy) in zip(ax2, SEASONS):
    for label, model, color in SATS:
        ax.plot(LON, eej(model, LON, [doy])[:, 0], color=color, lw=1.4, label=label)
    ax.axhline(0, color="0.7", lw=0.6)
    ax.text(0.02, 0.06, name, transform=ax.transAxes, fontsize=9,
            bbox=dict(fc="white", ec="none", alpha=0.7, pad=1.5))
    ax.set_ylabel("J (A/m)")
    ax.set_xlim(-180, 180)
    ax.set_xticks(range(-180, 181, 45))
    ax.set_ylim(-0.03, 0.18)
    ax.grid(True)
ax2[0].legend(fontsize=8, loc="upper right")
ax2[-1].set_xlabel(r"$\phi$ (degrees)")
fig2.tight_layout(rect=(0, 0, 1, 0.94))

# ----------------------------------------------------- Figure 5  (J vs season, 5 longitudes)
DOY = np.arange(0.0, 365.1, 3.0)
PHI = [-100, -45, 0, 45, 100]

fig5, ax5 = plt.subplots(len(PHI), 1, figsize=(5.5, 9), sharex=True)
fig5.suptitle("EEJ seasonal dependence by longitude\n(Alken & Maus 2007, Fig. 5)")
for ax, phi in zip(ax5, PHI):
    for label, model, color in SATS:
        ax.plot(DOY, eej(model, [phi], DOY)[0], color=color, lw=1.4, label=label)
    ax.axhline(0, color="0.7", lw=0.6)
    ax.text(0.97, 0.88, r"$\phi$ = %d°" % phi, transform=ax.transAxes, ha="right", va="top", fontsize=9)
    ax.set_ylabel("J (A/m)")
    ax.set_xlim(0, 365)
    ax.set_ylim(-0.01, 0.16)
    ax.grid(True)
ax5[0].legend(fontsize=8, loc="upper left")
ax5[-1].set_xlabel("day of year")
fig5.tight_layout(rect=(0, 0, 1, 0.94))

plt.show()
