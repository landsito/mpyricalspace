"""
Weimer ground geomagnetic perturbations -- Figure 8 of Weimer (2005b)
---------------------------------------------------------------------

"Geographic contour maps of predicted north (top) and east (bottom) components of the geomagnetic perturbations
at 2100 UT on 26 June 2001. The locations of six of the ground magnetometer stations are shown with the X marks.
The red line shows the lower-latitude boundary that is shared by the electric and magnetic potential models.
All modeled ionospheric currents flow inside this boundary, but the geomagnetic effects extend to lower
latitudes. The numbers in the lower left and right corners of the two maps show the minimum and maximum values."

Weimer, D. R. (2005), Predicting surface geomagnetic variations using ionospheric electrodynamic models,
J. Geophys. Res., 110, A12307, doi:10.1029/2005JA011270.

The model gives the electric potential Phi (`epot`) and the magnetic Euler potential psi/mu0 (`mpot`); the
ground perturbation is not part of it. This script follows the paper:

  1. E = -grad(Phi) and the potential current J_p = grad(psi)/mu0 (eq. 4) on an AACGM grid at 110 km.
  2. The equivalent (Hall) current from both potentials, eq. (10):  J_e = B^ x alpha (J_p . E^) E^,  with
     alpha = Sigma_H / Sigma_P = 1.5 as in the paper.
  3. Each grid cell's current element J_e dA is carried from AACGM to geographic coordinates with the local
     Jacobian of the mapping (the compensation for the non-uniform mapping the paper refers to, Richmond 1995).
  4. Biot-Savart from the current sheet at 110 km to the ground, eq. (12), using only the currents within an
     arc of 8 deg of the point overhead, as the paper does.

Drivers: the paper propagates ACE data to the magnetosphere with its own minimum-variance technique; here they
are NOAA OMNI (already shifted to the bow shock), averaged over the 20 min before 21:00 UT (the 5-min samples
of 20:45-21:00) the way 2005b drives the model. They are written below so that the script needs no index store.
The maps show the same quantities, at the same time and in the same layout as the paper's figure. This is an
approximate reproduction: the paper's own ACE propagation and coordinate mapping are not reproduced, so differences
from the published values are expected.

The paper's figure shows nothing equatorward of the model boundary (D. Weimer, pers. comm.: his graphics routines
took everything below the boundary as zero), although the ground perturbation can in principle be computed there.
MASK_OUTSIDE = True (the default) follows the paper and leaves those points blank; set it to False to draw the
perturbation at lower latitudes too.

    pip install "mpyricalspace[examples]"          # for matplotlib
    python examples/plot_weimer05_fig08.py
"""
import warnings
from datetime import datetime

import aacgmv2
import matplotlib.pyplot as plt
import numpy as np

from mpyricalspace import models

T = datetime(2001, 6, 26, 21)
# NOAA OMNI 5-min, mean of 20:45, 20:50, 20:55, 21:00 UT -- by [4.21, 1.57, -1.05, -1.28], bz [-5.99, -7.15, -7.11,
# -7.20] nT, vsw [473.7, 478.6, 483.7, 477.4] km/s, nsw [3.15, 3.19, 3.31, 3.47] cm^-3
DRV = dict(by=0.8625, bz=-6.8625, vsw=478.35, nsw=3.28)
ALPHA = 1.5                                         # Hall / Pedersen conductance ratio (2005b, sec. 4)
ARC_DEG = 8.0                                       # Biot-Savart only within this arc of the point overhead
RE, H = 6371.2e3, 110e3                             # m; the current sheet is at 110 km
MU0 = 4e-7 * np.pi
MASK_OUTSIDE = True                                 # blank the maps equatorward of the model boundary, as the paper

STATIONS = {"THL": (77.5, 290.8), "GDH": (69.2, 306.5), "NAQ": (61.2, 314.6), "FCC": (58.8, 265.9),
            "YKC": (62.5, 245.5), "BRW": (71.3, 203.4)}               # Table 1 of the paper (geographic)

# ------------------------------------------------------------------ 1. potentials on an AACGM grid at 110 km
D_DEG = 0.5
mlat = np.arange(40.0 + D_DEG / 2, 90.0, D_DEG)                       # cell centres (the pole is a cell edge)
mlon = np.arange(D_DEG / 2, 360.0, D_DEG)
LAT, LON = np.meshgrid(mlat, mlon, indexing="ij")
mlt = aacgmv2.convert_mlt(LON.ravel(), T, m2a=False)
ds = models.weimer05([T] * LAT.size, LAT.ravel(), mlt, tilt=None, **DRV)        # aligned samples
inside = np.isfinite(ds.epot.values).reshape(LAT.shape)
phi = np.nan_to_num(ds.epot.values.reshape(LAT.shape)) * 1e3                    # V; zero on/after the boundary
psi = np.nan_to_num(ds.mpot.values.reshape(LAT.shape)) * 1e3                    # A  (psi / mu0)


def grad(f):
    "(east, north) surface gradient on the 110 km sphere of a (mlat, mlon) field; periodic in longitude"
    r, h = RE + H, np.radians(D_DEG)
    north = np.gradient(f, h, axis=0) / r
    east = (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / (2 * h) / (r * np.cos(np.radians(LAT)))
    return east, north


# ------------------------------------------------------------------ 2. equivalent current, eq. (10)
ee, en = (-g for g in grad(phi))                                      # E = -grad(Phi)        [V/m]
pe, pn = grad(psi)                                                    # J_p = grad(psi)/mu0   [A/m], eq. (4)
emag = np.hypot(ee, en)
ue, un = (np.where(emag > 0, c / np.where(emag > 0, emag, 1), 0.0) for c in (ee, en))
proj = ALPHA * (pe * ue + pn * un)                                    # alpha (J_p . E^)
# B^ points down in the northern hemisphere: B^ x (v_e, v_n) = (v_n, -v_e)
je, jn = proj * un, -proj * ue                                        # J_e [A/m]

# ------------------------------------------------------------------ 3. current elements to geographic coordinates
DA = (RE + H) ** 2 * np.cos(np.radians(LAT)) * np.radians(D_DEG) ** 2    # cell area [m^2]
cell = inside & (np.hypot(je, jn) > 0)
la, lo = LAT[cell], LON[cell]


def a2g(lat, lon):
    glat, glon, _ = aacgmv2.convert_latlon_arr(lat, lon, H / 1e3, T, method_code="A2G")
    return glat, glon


def geo_step(g0, g1):
    "(east, north) distance [m] on the 110 km sphere between two geographic points (small steps)"
    dlon = (g1[1] - g0[1] + 180.0) % 360.0 - 180.0
    return ((RE + H) * np.radians(dlon) * np.cos(np.radians(g0[0])), (RE + H) * np.radians(g1[0] - g0[0]))


DELTA = 0.05                                                          # deg, finite-difference step of the Jacobian
g0 = a2g(la, lo)
dn_e, dn_n = geo_step(g0, a2g(la + DELTA, lo))                        # a magnetic step north ...
de_e, de_n = geo_step(g0, a2g(la, lo + DELTA / np.cos(np.radians(la))))   # ... and one east, of equal length
step = (RE + H) * np.radians(DELTA)
mm_e, mm_n = je[cell] * DA[cell], jn[cell] * DA[cell]                 # current elements [A m], magnetic frame
m_e = (de_e * mm_e + dn_e * mm_n) / step                              # the Jacobian carries them to geographic
m_n = (de_n * mm_e + dn_n * mm_n) / step
ok = np.isfinite(g0[0]) & np.isfinite(m_e) & np.isfinite(m_n)
glat, glon, m_e, m_n = g0[0][ok], g0[1][ok], m_e[ok], m_n[ok]


def unit_vectors(lat, lon):
    "radial, east and north unit vectors (Earth-centred cartesian) at geographic lat, lon [deg]"
    la_, lo_ = np.radians(lat), np.radians(lon)
    up = np.stack([np.cos(la_) * np.cos(lo_), np.cos(la_) * np.sin(lo_), np.sin(la_)], -1)
    east = np.stack([-np.sin(lo_), np.cos(lo_), np.zeros_like(lo_)], -1)
    north = np.stack([-np.sin(la_) * np.cos(lo_), -np.sin(la_) * np.sin(lo_), np.cos(la_)], -1)
    return up, east, north


s_up, s_e, s_n = unit_vectors(glat, glon)
r1 = (RE + H) * s_up                                                  # source positions
m = m_e[:, None] * s_e + m_n[:, None] * s_n                           # current elements, cartesian [A m]


# ------------------------------------------------------------------ 4. Biot-Savart to the ground, eq. (12)
def delta_b(lat, lon, chunk=200):
    "ground perturbation (north, east, down) [nT] at geographic points"
    o_up, o_e, o_n = unit_vectors(np.ravel(lat), np.ravel(lon))
    out = np.zeros((o_up.shape[0], 3))
    cos_arc = np.cos(np.radians(ARC_DEG))
    for i in range(0, o_up.shape[0], chunk):
        u = o_up[i:i + chunk]
        near = (u @ s_up.T) >= cos_arc                                # sources within ARC_DEG of the point overhead
        d = RE * u[:, None, :] - r1[None, :, :]                       # r2 - r1
        b = np.cross(m[None, :, :], d) / np.linalg.norm(d, axis=-1, keepdims=True) ** 3
        b = MU0 / (4 * np.pi) * np.where(near[..., None], b, 0.0).sum(axis=1) * 1e9      # nT
        out[i:i + chunk] = np.stack([(b * o_n[i:i + chunk]).sum(-1), (b * o_e[i:i + chunk]).sum(-1),
                                     -(b * u).sum(-1)], -1)
    return out.reshape(np.shape(lat) + (3,))


glat_map = np.arange(40.0, 90.01, 1.0)
glon_map = np.arange(0.0, 360.01, 3.0)
GLA, GLO = np.meshgrid(glat_map, glon_map, indexing="ij")
db = delta_b(GLA, GLO)

# the model boundary (the lowest latitude inside the cap, for each magnetic longitude), in geographic
first = np.argmax(inside, axis=0)
blat, blon = a2g(mlat[first] - D_DEG / 2, mlon)
if MASK_OUTSIDE:
    # a ground point is inside when the point overhead (110 km) lies poleward of the boundary in AACGM
    pla, plo, _ = aacgmv2.convert_latlon_arr(GLA.ravel(), GLO.ravel(), H / 1e3, T, method_code="G2A")
    edge = np.interp(np.mod(plo, 360.0), np.r_[mlon - 360.0, mlon, mlon + 360.0], np.tile(mlat[first] - D_DEG / 2, 3))
    outside = ~(np.isfinite(pla) & (pla >= edge)).reshape(GLA.shape)
    db[outside] = np.nan

# ------------------------------------------------------------------ plot: polar maps, as oriented in the paper
# The paper's maps: an orthographic view of the north pole with the 225 deg E meridian at the top (NAQ to the left,
# BRW above to the right), parallels every 10 deg and meridians every 15 deg. With cartopy installed the maps carry
# the continents and country borders (Natural Earth, downloaded by cartopy on first use); without it, plain polar
# axes (no map).
TOP_LON = 225.0
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except ImportError:
    ccrs = None
    warnings.warn("plot_weimer05_fig08: cartopy is not installed, so the maps are drawn without continents or "
                  "borders (pip install cartopy to add them)")


def map_axes(fig, pos):
    "one polar map; returns the axes and a function (lon, lat) -> the keyword arguments to plot there"
    if ccrs is not None:
        ax = fig.add_subplot(pos, projection=ccrs.Orthographic(central_longitude=TOP_LON - 180.0, central_latitude=90))
        r40 = ax.projection.globe.semimajor_axis or 6378137.0
        r40 *= np.cos(np.radians(40.0))                               # the 40 deg parallel, in projected metres
        ax.set_xlim(-r40, r40)                                        # (set_extent does not centre a polar view)
        ax.set_ylim(-r40, r40)
        theta = np.linspace(0, 2 * np.pi, 200)
        ax.set_boundary(plt.matplotlib.path.Path(np.c_[np.sin(theta), np.cos(theta)] * 0.5 + 0.5), transform=ax.transAxes)
        ax.add_feature(cfeature.LAND.with_scale("110m"), facecolor="0.85", zorder=0)
        ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.6, zorder=3)
        ax.add_feature(cfeature.BORDERS.with_scale("110m"), linewidth=0.4, linestyle=":", zorder=3)
        ax.gridlines(ylocs=[40, 50, 60, 70, 80], xlocs=np.arange(-180, 180, 15), linewidth=0.4, color="0.5")
        return ax, lambda lon, lat: (np.asarray(lon), np.asarray(lat), {"transform": ccrs.PlateCarree()})
    ax = fig.add_subplot(pos, projection="polar")
    ax.set_theta_offset(np.radians(90.0 - TOP_LON))
    ax.set_ylim(0, 50)
    ax.set_yticks([10, 20, 30, 40])
    ax.set_yticklabels(["80", "70", "60", "50"], fontsize=7)
    ax.set_xticks(np.radians(np.arange(0, 360, 15)))
    ax.set_xticklabels([])
    return ax, lambda lon, lat: (np.radians(lon), 90.0 - np.asarray(lat), {})


fig = plt.figure(figsize=(6.5, 11))
fig.suptitle("Weimer ground $\\Delta$B, %s UT (Weimer 2005b, Fig. 8)\nBy %.1f, Bz %.1f nT, V %.0f km/s, N %.1f cm$^{-3}$"
             % (T.strftime("%d %b %Y %H:%M"), DRV["by"], DRV["bz"], DRV["vsw"], DRV["nsw"]))
for pos, (name, k, step_nt, top) in zip((211, 212), [("N-S", 0, 35, 280), ("E-W", 1, 25, 200)]):
    ax, at = map_axes(fig, pos)
    z = db[..., k]
    levels = np.arange(-top, top + 0.1, step_nt)
    x, y, kw = at(GLO, GLA)
    cf = ax.contourf(x, y, z, levels=levels, cmap="RdBu_r", extend="both", alpha=0.8, zorder=1, **kw)
    ax.contour(x, y, z, levels=levels, colors="k", linewidths=0.3, zorder=2, **kw)
    x, y, kw = at(np.append(blon, blon[0]), np.append(blat, blat[0]))
    if kw:
        kw = {"transform": ccrs.Geodetic()}                          # a line across the date line stays on the map
    ax.plot(x, y, color="r", lw=1.2, zorder=4, **kw)
    for code, (sla, slo) in STATIONS.items():
        x, y, kw = at(slo, sla)
        ax.plot(x, y, "kx", ms=6, zorder=5, **kw)
        ax.text(x, y, "  " + code, fontsize=7, va="bottom", zorder=5, **kw)
    lo_, hi_ = float(np.nanmin(z)), float(np.nanmax(z))
    ax.annotate("%.0f nT" % lo_, (0.0, -0.04), xycoords="axes fraction", ha="left", fontsize=9)
    ax.annotate("%.0f nT" % hi_, (1.0, -0.04), xycoords="axes fraction", ha="right", fontsize=9)
    fig.colorbar(cf, ax=ax, label="%s $\\Delta$B [nT]" % name, shrink=0.8)
plt.show()
