"""
Weimer high-latitude electric potential -- Figure 2 of Weimer (2005a)
---------------------------------------------------------------------

"Polar cap electric potentials in the Northern Hemisphere, mapped as a function of AACGM latitude and MLT.
Panels a-d and f-i show the patterns for eight different clock angle orientations of the IMF vector in the
GSM Y-Z plane; the angle in degrees is indicated in the top left corner of each map. The IMF has a fixed
magnitude of 5 nT, the solar wind velocity is 450 km/s, the solar wind number density is 4 cm^-3, and the
dipole tilt angle is 0. Panel e shows the potential for zero IMF, with the same solar wind conditions.
Minimum and maximum potential values are printed in the bottom left and right corners of each map."

Weimer, D. R. (2005), Improved ionospheric electrodynamic models and application to calculating Joule
heating rates, J. Geophys. Res., 110, A05306, doi:10.1029/2004JA010884.

The paper shows the earlier hybrid version of the models; the spherical-cap-harmonic revision that
mpyricalspace bundles reproduces its potentials closely. Drivers are given explicitly here, so nothing is
pulled from the index store. The script checks the extremes against the values printed in the paper's
figure and fails if they drift by more than TOL_KV.

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_weimer05_fig02.py
"""
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt

from mpyricalspace import Predictor

BT, VSW, NSW, TILT = 5.0, 450.0, 4.0, 0.0
TOL_KV = 4.0                                            # the paper is the older hybrid model: a benchmark, not a target

# (panel, IMF clock angle [deg] -- None = zero IMF, min/max potential printed in the paper [kV])
PANELS = [("a", 315, (-17, 15)), ("b", 0, (-12, 12)), ("c", 45, (-22, 18)),
          ("d", 270, (-30, 28)), ("e", None, (-18, 15)), ("f", 90, (-40, 29)),
          ("g", 225, (-46, 44)), ("h", 180, (-55, 47)), ("i", 135, (-52, 41))]

MLAT = np.arange(45.0, 90.01, 0.5)
MLT = np.arange(0.0, 24.0, 0.25)
LEVELS = np.arange(-54.0, 54.1, 6.0)

# one call, nine "times" (the drivers differ per panel): a (time, mlat, mlt) cube
theta = np.radians([0.0 if c is None else c for _, c, _ in PANELS])
bt = np.array([0.0 if c is None else BT for _, c, _ in PANELS])
m = Predictor.Empirical()
m.set_time([datetime(2024, 3, 20, 12) + timedelta(minutes=i) for i in range(len(PANELS))])
ds = m.run_weimer05(MLAT, MLT, by=bt * np.sin(theta), bz=bt * np.cos(theta), vsw=VSW, nsw=NSW, tilt=TILT)

# polar axes as in the paper: 12 MLT at the top, dawn (6 MLT) on the right, radius = 90 - AACGM latitude
# (zero angle at the bottom = 0 MLT, angle increasing counter-clockwise = increasing MLT)
ang = np.radians(np.append(MLT, 24.0) * 15.0)
TH, R = np.meshgrid(ang, 90.0 - MLAT)

fig, axes = plt.subplots(3, 3, figsize=(10, 9.5), subplot_kw={"projection": "polar"})
fig.suptitle("Weimer electric potential, BT = %.0f nT, V = %.0f km/s, N = %.0f cm$^{-3}$, tilt %.0f$^\\circ$\n"
             "(Weimer 2005, Fig. 2)" % (BT, VSW, NSW, TILT))
cf, worst = None, 0.0
print("panel  clock   min/max [kV] paper -> here")
for ax, (name, clock, (pmin, pmax)), pot in zip(axes.ravel(), PANELS, ds.epot.values):
    z = np.nan_to_num(pot, nan=0.0)                     # zero on and beyond the model's low-latitude boundary
    z = np.hstack([z, z[:, :1]])                        # close the circle at 24 MLT
    cf = ax.contourf(TH, R, z, levels=LEVELS, cmap="RdYlBu_r", extend="both")
    ax.contour(TH, R, z, levels=LEVELS, colors="k", linewidths=0.3)
    ax.set_theta_zero_location("S")
    ax.set_ylim(0, 45)
    ax.set_rlabel_position(157.5)
    ax.set_yticks([10, 20, 30, 40])
    ax.set_yticklabels(["80", "70", "60", "50"], fontsize=6)
    ax.set_xticks(np.radians(np.arange(0, 24, 2) * 15.0))
    ax.set_xticklabels([str(h) for h in range(0, 24, 2)], fontsize=7)
    lo, hi = float(np.nanmin(pot)), float(np.nanmax(pot))
    ax.set_title("%s  %s" % (name, "zero IMF" if clock is None else "%d$^\\circ$" % clock), loc="left", fontsize=9)
    ax.annotate("%.0f kV" % lo, (0.0, -0.06), xycoords="axes fraction", ha="left", fontsize=8)
    ax.annotate("%.0f kV" % hi, (1.0, -0.06), xycoords="axes fraction", ha="right", fontsize=8)
    print("  %s   %-8s %4d %4d  ->  %6.1f %5.1f" % (name, "zero" if clock is None else clock, pmin, pmax, lo, hi))
    worst = max(worst, abs(lo - pmin), abs(hi - pmax))

fig.subplots_adjust(right=0.88, hspace=0.3, wspace=0.3, top=0.9)
cax = fig.add_axes([0.91, 0.25, 0.02, 0.5])
fig.colorbar(cf, cax=cax, label="kV")
print("largest difference from the paper's printed extremes: %.1f kV (allowed %.1f)" % (worst, TOL_KV))
assert worst <= TOL_KV, "Weimer05 potentials drifted from Weimer (2005a) Fig. 2 by %.1f kV" % worst
plt.show()
