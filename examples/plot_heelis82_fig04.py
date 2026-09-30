"""
Heelis convection potential -- Figure 4 of Heelis, Lowell and Spiro (1982), and the TIE-GCM version
-----------------------------------------------------------------------------------------------------

"Figure 4 shows the local time distribution of potential along the convection reversal boundary, assuming that
the dawnside boundary potential is +20 kV, the duskside boundary potential is -35 kV, phi_d = 15, phi_n = 180,
phi_d+ = phi_d- = 22.5, phi_n+ = phi_n- = 55. [...] Regions of electrostatic equipotential exist between 13:00 hours
and 21:00 hours and between 04:00 hours and 09:00 hours when there is no flow across the boundary."

Heelis, R. A., J. K. Lowell, and R. W. Spiro (1982), A model of the high-latitude ionospheric convection pattern,
J. Geophys. Res., 87(A8), 6339-6345, doi:10.1029/JA087iA08p06339.

mpyricalspace's `heelis82` model is NCAR TIE-GCM 2.0's implementation of this pattern (the default `variant='tiegcm'`
reproduces what TIE-GCM computes with its Heelis potential model). `variant='paper'` swaps only TIE-GCM's latitude
function for the paper's, and with `params=` given the paper's own constants it reproduces Fig. 4. The four panels:

  a  Fig. 4 of the paper: potential along the convection reversal boundary (variant='paper', the paper's constants).
     The script checks the plateaus (+20 and -35 kV, at the local times the paper gives) and the potential at the
     centres of the two convergence zones (the mean, -7.5 kV) and fails if any is off by more than TOL_KV.
  b  Potential across the boundary on the dusk meridian, for both variants: TIE-GCM's has a kink (the north-south
     electric field jumps), the paper's elliptical smoothing removes it. The script checks that.
  c, d  The pattern as TIE-GCM computes it (c) and with the paper's latitude function (d), for Kp = 3 (the cross-cap
     potential follows from Kp: 15 + 15 Kp + 0.8 Kp^2 = 67.2 kV) and By = 0. The script checks that the potential
     drop across the cap (max - min) is that cross-cap potential.

Kp, By and the constants are given explicitly here, so nothing is pulled from the index store.

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_heelis82_fig04.py
"""
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from mpyricalspace import Predictor

TOL_KV = 0.1                                           # the check against the paper's stated values
KP = 3.0                                               # -> cross-cap potential 67.2 kV
T = datetime(2024, 5, 11, 12)                          # the drivers are given: the date only labels the result

# Fig. 4 constants. phi_d and phi_n are hour angles clockwise from noon in the paper; here local times
# (phi_d = 15 deg -> 11 MLT, phi_n = 180 deg -> 24 MLT). theta0 = 15 deg; the pattern is centred on the pole.
PAPER = dict(theta0=15., psim=20., psie=-35., pcen=0., phid=11., phin=24., phidp=22.5, phidm=22.5, phinp=55., phinm=55.,
             offc=0., dskofc=0.)
# (local time [h], potential [kV]) the paper's figure and text give
PLATEAUS = [(4.0, 20.), (7.0, 20.), (9.0, 20.), (13.0, -35.), (17.0, -35.), (20.0, -35.)]   # 04-09 h and 13-21 h
CENTRES = [(11.0, -7.5), (0.0, -7.5)]                  # phi_d and phi_n: half way between +20 and -35

m = Predictor.Empirical()
m.set_time([T])

# ---- a: Fig. 4, along the reversal boundary (colatitude theta0 = 15 deg, just outside so the outer branch is used)
MLT = np.arange(0.0, 24.0, 0.05)
bnd = m.run_heelis82([90. - 15.0001], MLT, cp=45., by=0., variant="paper", params=PAPER).epot.values

# ---- b: across the boundary on the dusk meridian (18 MLT); pattern centred on the pole so colatitude = 90 - latitude
LAT = np.arange(70., 82.001, 0.05)
CENTRED = dict(offc=0., dskofc=0.)
prof = {v: m.run_heelis82(LAT, [18.], cp=45., by=0., variant=v, params=CENTRED).epot.values.ravel()
        for v in ("tiegcm", "paper")}
jump = {v: np.abs(np.diff(np.diff(p) / 0.05)).max() for v, p in prof.items()}    # largest change of dV/dlat [kV/deg]

# ---- c, d: polar maps
MLAT = np.arange(50.0, 90.01, 0.1)
MAP_MLT = np.arange(0.0, 24.0, 0.1)
maps = {v: m.run_heelis82(MLAT, MAP_MLT, kp=KP, by=0., variant=v) for v in ("tiegcm", "paper")}

# ---- checks
print("panel a   local time  paper -> here [kV]")
worst = 0.0
for hh, want in PLATEAUS + CENTRES:
    got = float(bnd[int(round(hh / 0.05)) % MLT.size])
    print("          %5.1f h  %7.1f -> %7.2f" % (hh, want, got))
    worst = max(worst, abs(got - want))
print("panel a   largest difference from the paper's values: %.3f kV (allowed %.1f)" % (worst, TOL_KV))
print("panel b   largest change of dV/dlat between adjacent 0.05 deg cells: TIE-GCM %.2f, paper %.2f kV/deg"
      % (jump["tiegcm"], jump["paper"]))
cp_kp = 15. + 15. * KP + 0.8 * KP ** 2
for v, ds in maps.items():
    drop = float(ds.epot.max() - ds.epot.min())
    print("panel %s   %-7s cross-cap potential: max - min = %.2f kV (cp from Kp = %.2f kV, ds.cp = %.2f)"
          % ("c" if v == "tiegcm" else "d", v, drop, cp_kp, float(ds.cp)))
    assert abs(drop - cp_kp) < 0.5, "%s: the potential drop across the cap is not the cross-cap potential" % v
assert worst <= TOL_KV, "the potential along the boundary differs from Heelis et al. (1982) Fig. 4 by %.3f kV" % worst
assert jump["tiegcm"] > 5 * jump["paper"], "the TIE-GCM kink / the paper's smoothing were not found"

# ---- plot
fig = plt.figure(figsize=(11, 9))
ax = fig.add_subplot(2, 2, 1)
ax.plot(MLT, bnd, "k")
for a, b in ((4, 9), (13, 21)):                                            # the equipotential regions of the paper
    ax.axvspan(a, b, color="0.9", zorder=0)
for hh, want in CENTRES:
    ax.plot([hh], [want], "o", ms=4, color="C3")
ax.set(xlim=(0, 24), xticks=range(0, 25, 3), xlabel="MLT [h]", ylabel="potential [kV]",
       title="a  Fig. 4: along the reversal boundary")
ax.text(0.02, 0.05, "grey: equipotential regions the paper gives;\nred: centres of the convergence zones",
        transform=ax.transAxes, fontsize=7)

ax = fig.add_subplot(2, 2, 2)
for v, c in (("tiegcm", "C3"), ("paper", "C0")):
    ax.plot(LAT, prof[v], c, label="%s  (max dV/dlat change %.2f kV/deg)" % (v, jump[v]))
ax.set(xlabel="latitude [deg], 18 MLT", ylabel="potential [kV]", title="b  across the boundary (cross-cap 45 kV)")
ax.legend(fontsize=7, loc="lower right")

LEVELS = np.arange(-40.0, 35.1, 5.0)
ang = np.radians(np.append(MAP_MLT, 24.0) * 15.0)
TH, R = np.meshgrid(ang, 90.0 - MLAT)
cf = None
for k, (v, title) in enumerate((("tiegcm", "c  TIE-GCM (variant='tiegcm')"),
                                ("paper", "d  paper's latitude function (variant='paper')"))):
    ax = fig.add_subplot(2, 2, 3 + k, projection="polar")
    z = maps[v].epot.values
    z = np.hstack([z, z[:, :1]])                                           # close the circle at 24 MLT
    cf = ax.contourf(TH, R, z, levels=LEVELS, cmap="RdYlBu_r", extend="both")
    ax.contour(TH, R, z, levels=LEVELS, colors="k", linewidths=0.3)
    ax.set_theta_zero_location("S")                                         # 12 MLT at the top, dawn (6 MLT) on the right
    ax.set_ylim(0, 40)
    ax.set_rlabel_position(157.5)
    ax.set_yticks([10, 20, 30, 40])
    ax.set_yticklabels(["80", "70", "60", "50"], fontsize=6)
    ax.set_xticks(np.radians(np.arange(0, 24, 2) * 15.0))
    ax.set_xticklabels([str(hh) for hh in range(0, 24, 2)], fontsize=7)
    ax.set_title(title, loc="left", fontsize=9)
    ax.annotate("%.0f kV" % z.min(), (0.0, -0.06), xycoords="axes fraction", ha="left", fontsize=8)
    ax.annotate("%.0f kV" % z.max(), (1.0, -0.06), xycoords="axes fraction", ha="right", fontsize=8)
fig.suptitle("Heelis convection potential (Heelis et al. 1982; TIE-GCM 2.0 implementation), Kp = %.0f, By = 0" % KP)
fig.subplots_adjust(right=0.9, hspace=0.35, wspace=0.3, top=0.92)
cax = fig.add_axes([0.93, 0.1, 0.015, 0.3])
fig.colorbar(cf, cax=cax, label="kV")
plt.show()
