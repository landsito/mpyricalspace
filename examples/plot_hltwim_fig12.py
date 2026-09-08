"""
HL-TWiM high-latitude thermospheric winds -- Figure 12 of Dhadly et al. (2019)
-----------------------------------------------------------------------------

"Polar wind vector plots for various UTs (0, 8 and 16) and Kp (2 and 5) for
northern (top two rows) and southern (bottom two rows) high latitudes as a
function of MLAT and MLT."

A regular magnetic-latitude / magnetic-local-time grid is inverted to geographic
coordinates (apexpy) for each UT, HL-TWiM is run there, and the quasi-dipole
winds are drawn on an MLT polar dial.

    pip install "mpyricalspace[examples]"          # matplotlib + apexpy
    python examples/plot_hltwim_fig12.py [--doy 80] [--scale 6]
"""
import argparse
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt

try:
    import apexpy
except ImportError:
    raise SystemExit('this example needs apexpy:  pip install "mpyricalspace[examples]"')

from mpyricalspace import Predictor

ap = argparse.ArgumentParser()
ap.add_argument("--doy", type=int, default=200, help="day of year (default 200)")
ap.add_argument("--year", type=int, default=2016)
ap.add_argument("--scale", type=float, default=18.0, help="m/s per plot degree (smaller -> longer arrows)")
args = ap.parse_args()

UTS = (0, 8, 16)
MLAT = np.arange(41.0, 88.1, 5)                 # magnetic-latitude rings
MLT = np.arange(0.0, 24.0, 1.0)                   # magnetic-local-time sectors
GLA, GLT = np.meshgrid(MLAT, MLT, indexing="ij")
R = 90.0 - GLA
ANG = np.deg2rad(GLT * 15.0 - 90.0)              # 00 MLT bottom, 06 right, 12 top, 18 left
RMAX = 90.0 - MLAT[0]


def winds(hemi, kp, ut):
    """QD winds (mu east, mv north) on the MLAT/MLT grid for one hemisphere / Kp / UT."""
    dt = datetime(args.year, 1, 1) + timedelta(days=args.doy - 1, hours=ut)
    A = apexpy.Apex(date=dt)
    mlon = A.mlt2mlon(GLT.ravel(), dt)
    glat, glon, _ = A.qd2geo(hemi * GLA.ravel(), mlon, 250.0)
    m = Predictor.Empirical()
    m.set_time([dt])
    ds = m.run_hltwim(kp=float(kp), lat=glat, lon=glon, ut=float(ut), doy=args.doy)
    return ds.mu.data.reshape(GLA.shape), ds.mv.data.reshape(GLA.shape)


def panel(ax, hemi, kp, ut):
    mu, mv = winds(hemi, kp, ut)
    px, py = R * np.cos(ANG), R * np.sin(ANG)
    v_out = -mv if hemi > 0 else mv                          # QD north: inward (N) / outward (S)
    dx = mu * (-np.sin(ANG)) + v_out * np.cos(ANG)
    dy = mu * (np.cos(ANG)) + v_out * np.sin(ANG)

    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-RMAX - 6, RMAX + 6); ax.set_ylim(-RMAX - 6, RMAX + 6)
    th = np.linspace(0, 2 * np.pi, 180)
    for lat in (40, 50, 60, 70, 80):
        rr = 90.0 - lat
        ax.plot(rr * np.cos(th), rr * np.sin(th), color="0.6",
                lw=(1.0 if lat == 40 else 0.4), ls=("solid" if lat == 40 else (0, (4, 3))))
        ax.text(rr, -1.4, str(lat), fontsize=5.5, color="0.45", ha="center", va="top")
    for a in np.deg2rad([-90, 0, 90, 180]):
        ax.plot([0, RMAX * np.cos(a)], [0, RMAX * np.sin(a)], color="0.6", lw=0.4, ls=(0, (4, 3)))
    ax.plot(0, 0, "+", color="0.35", ms=6)
    e = RMAX + 3
    for txt, (x, y, ha, va) in {"12": (0, e, "center", "bottom"), "00": (0, -e, "center", "top"),
                                "06": (e, 0, "left", "center"), "18": (-e, 0, "right", "center")}.items():
        ax.text(x, y, txt, ha=ha, va=va, fontsize=7)

    q = ax.quiver(px, py, dx, dy, color="tab:blue", angles="xy", scale_units="xy",
                  scale=args.scale, width=0.005, headwidth=3.5, headlength=4)
    ax.set_title("UT=%d, KP=%d" % (ut, kp), fontsize=9)
    return q


def main():
    fig, axes = plt.subplots(4, 3, figsize=(7, 8.8))
    fig.suptitle("HL-TWiM polar winds\n(Dhadly et al. 2019, Fig. 12;  doy = %d)" % args.doy)
    rows = [(1, 2), (1, 5), (-1, 2), (-1, 5)]               # (hemisphere, Kp) per row
    q = None
    for ri, (hemi, kp) in enumerate(rows):
        for ci, ut in enumerate(UTS):
            q = panel(axes[ri, ci], hemi, kp, ut)

    fig.text(0.035, 0.74, "Northern hemisphere", rotation=90, va="center", fontsize=10)
    fig.text(0.035, 0.30, "Southern hemisphere", rotation=90, va="center", fontsize=10)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.84, bottom=0.02, hspace=0.35, wspace=0.05)
    axes[0, 1].quiverkey(q, 0.5, 0.9, 100, "100 m/s", labelpos="E",
                         coordinates="figure", color="tab:blue")
    plt.show()


if __name__ == "__main__":
    main()
