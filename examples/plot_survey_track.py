"""
Many models at once along a low-inclination satellite track
----------------------------------------------------------

``mpyricalspace.survey.run_track`` evaluates every applicable empirical model at a
set of (time, lat, lon, alt) samples and merges the result into one xr.Dataset.
The global models (IGRF, MSIS, HWM, IRI) run at every point; the
equatorial-electrodynamics models (``sf``, ``eef``, ``manoj`` ...) appear only
where the track is within a few degrees of the dip equator and inside their
altitude regime -- so on a low-inclination orbit they light up twice per orbit,
at the northbound / southbound equator crossings.

The ground track here is analytic (no ephemeris file): for a circular orbit of
inclination ``i`` and period ``T``, with ``u = 2*pi*t/T`` the argument of latitude,

    lat(t) = asin( sin(i) * sin(u) )                       -- a latitude sinusoid
    lon(t) = lon0 + atan2(cos(i)*sin(u), cos(u)) - 360*t/1436    -- Earth turns under the orbit

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_survey_track.py
"""
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt

from mpyricalspace import survey

# --- an analytic low-inclination orbit -------------------------------------
EPOCH   = datetime(2004, 11, 8, 0, 0)                   # 2004 -> local RTEEF (manoj) works
INC     = 20.0                                          # inclination [deg]
ALT     = 480.0                                         # km  (circular)
T_ORBIT = 94.4                                          # min (~480 km)
N_ORBIT = 3
DT      = 1.0                                           # min sampling
LON0    = -60.0                                         # sub-satellite longitude at the epoch

tmin = np.arange(0.0, N_ORBIT * T_ORBIT, DT)
u    = 2 * np.pi * tmin / T_ORBIT
lat  = np.degrees(np.arcsin(np.sin(np.radians(INC)) * np.sin(u)))
lon  = np.degrees(np.arctan2(np.cos(np.radians(INC)) * np.sin(u), np.cos(u)))
lon  = (LON0 + lon - 360.0 * tmin / 1436.0 + 180.0) % 360.0 - 180.0        # + Earth rotation
alt  = np.full(tmin.size, ALT)
times = [EPOCH + timedelta(minutes=float(m)) for m in tmin]
th   = tmin / 60.0                                      # hours, for the x axis

# --- one call -------------------------------------------------------------
ds = survey.run_track(times, lat, lon, alt)
print("ran    :", ds.attrs["models"])
print("skipped:", ds.attrs["skipped"] or "-")

near_eq = np.abs(ds.dip_lat.data) <= 2.5                # the tightest equatorial band (sf)

fig, ax = plt.subplots(4, 1, figsize=(7, 9.5), constrained_layout=True)
fig.suptitle("Model state along a low-inclination orbit  (i=%g deg, %g km, %s)\n"
             "one survey.run_track() call" % (INC, ALT, EPOCH.date()))

# 1. ground track, coloured by time; equator crossings marked
sc = ax[0].scatter(lon, lat, c=th, cmap="viridis", s=6)
ax[0].plot(lon[near_eq], lat[near_eq], "r.", ms=4, label="|dip lat| <= 2.5 deg")
ax[0].axhline(0, color="0.7", lw=0.6)
ax[0].set_xlabel("longitude (deg)"); ax[0].set_ylabel("latitude (deg)")
ax[0].set_xlim(-180, 180); ax[0].set_xticks(range(-180, 181, 60))
ax[0].legend(fontsize=7, loc="upper right")
fig.colorbar(sc, ax=ax[0], label="time (h)", pad=0.01)

# 2. IRI along the track
ax[1].plot(th, ds.iri_NMF2.data / 1e11, "k-", lw=1.2)
ax[1].set_ylabel(r"NmF2 ($10^{11}\,m^{-3}$)")
a1 = ax[1].twinx()
a1.plot(th, ds.iri_HMF2.data, color="tab:orange", lw=1.0)
a1.set_ylabel("hmF2 (km)", color="tab:orange")
ax[1].set_title("IRI-2020  peak density / height", fontsize=9)

# 3. MSIS density + HWM wind speed
ax[2].semilogy(th, ds.msis_rho.data, "k-", lw=1.2)
ax[2].set_ylabel(r"MSIS $\rho$ (kg m$^{-3}$)")
a2 = ax[2].twinx()
a2.plot(th, np.hypot(ds.hwm_u.data, ds.hwm_v.data), color="tab:green", lw=1.0)
a2.set_ylabel("HWM |wind| (m/s)", color="tab:green")
ax[2].set_title("MSIS 2.1 density  /  HWM14 wind speed", fontsize=9)

# 4. equatorial models -- finite only near the dip-equator crossings
ax[3].axhline(0, color="0.7", lw=0.6)
l1, = ax[3].plot(th, ds.sf_qvdrift.data, "o-", color="tab:blue", ms=4, lw=0.8,
                 label="sf  quiet drift (m/s)")
ax[3].set_ylabel("vertical drift (m/s)", color="tab:blue")
a3 = ax[3].twinx()
l2, = a3.plot(th, ds.manoj_ppef.data, "s-", color="tab:red", ms=4, lw=0.8,
              label="manoj  RTEEF prompt E (mV/m)")
l3, = a3.plot(th, ds.eef.data, "^-", color="tab:purple", ms=4, lw=0.8,
              label="eef  quiet climatology (mV/m)")
a3.set_ylabel("equatorial zonal E (mV/m)", color="tab:red")
ax[3].set_title("equatorial electrodynamics -- only near the dip-equator crossings\n"
                "(Nov 2004 superstorm: RTEEF prompt penetration >> quiet climatology)", fontsize=9)
ax[3].legend(handles=[l1, l2, l3], fontsize=7, loc="upper left")

for a in ax[1:]:
    a.set_xlim(0, th[-1]); a.grid(True, axis="x", alpha=0.3)
ax[3].set_xlabel("time (h)")
plt.show()
