# Heelis high-latitude convection potential, as computed by TIE-GCM 2.0

> **What this is for.** The purpose of this folder is to **reproduce what NCAR's TIE-GCM 2.0 computes when it is run
> with its Heelis potential model** (`potential_model = 'HEELIS'`): the same source file (`heelis.F`), the same
> parameters, the same dependence on the cross-cap potential and IMF By, the same Kp relation. It is **not** a fresh
> implementation of the paper. The functional form is that of Heelis, Lowell and Spiro (1982); everything that makes
> it a model of the real ionosphere (the constants, and how they follow ctpoten and By) is TIE-GCM's tuning. An
> optional `variant='paper'` runs the paper's own latitude function instead (see [Two variants](#two-variants)).

> Heelis, R. A., J. K. Lowell, and R. W. Spiro (1982), A model of the high-latitude ionospheric convection pattern,
> J. Geophys. Res., 87(A8), 6339-6345, [doi:10.1029/JA087iA08p06339](https://doi.org/10.1029/JA087iA08p06339)
>
> Qian, L., et al. (2014), The NCAR TIE-GCM: A community model of the coupled thermosphere/ionosphere system, in
> *Modeling the Ionosphere-Thermosphere System*, AGU Geophysical Monograph 201, 73-83,
> [doi:10.1002/9781118704417.ch7](https://doi.org/10.1002/9781118704417.ch7)

The model is `models.heelis82` (and `Empirical.run_heelis82`); the compiled extension is `heelis82f`. The name
is the year of the paper.

## Use

```python
import numpy as np
from datetime import datetime
from mpyricalspace import models

mlat, mlt = np.arange(50., 90.1, 1.), np.arange(0., 24., 1.)          # magnetic latitude [deg], local time [h]
t = datetime(2024, 5, 11, 12)

ds = models.heelis82(t, mlat, mlt)                    # Kp and By from the DataManager (a (mlat, mlt) cube)
ds = models.heelis82(t, mlat, mlt, cp=60.)            # cross-cap potential given: Kp is not looked up
ds = models.heelis82(t, mlat, mlt, kp=5., by=-3.)     # or Kp and By by hand
ds.epot                                               # kV; + in the morning cell, - in the evening cell
ds.cp, ds.kp, ds.by                                   # the drivers used, along time
ds = obj.run_heelis82(mlat, mlt)                      # the same through the Empirical facade (obj.time)
```

| Input | Meaning |
|---|---|
| `mlat`, `mlt` | magnetic latitude [deg] (negative = southern hemisphere) and local time [h]. The natural choice is AACGM (`aacgmv2`); see [Differences from a TIE-GCM run](#differences-from-a-tie-gcm-run). If time, mlat and mlt have the same length they are aligned samples (a track), otherwise the result is a (time, mlat, mlt) cube. |
| `cp` | cross-cap potential [kV]. **Optional.** If given, it wins over Kp. |
| `kp` | Kp index (0-9), used only where `cp` is not given. `None` -> the DataManager (GFZ 3-hourly), interpolated in time as TIE-GCM does. |
| `by` | IMF By [nT] (GSM). `None` -> the DataManager (NOAA OMNI). Limited to -11..+7 nT, as TIE-GCM does. |
| `res`, `avg`, `lag` | how `by` is taken from the store, as in `weimer05`; default `res='5min', avg=20` (below). |
| `kp_interp` | `True` (default): Kp interpolated between the centres of its 3-h intervals, like TIE-GCM. `False`: the value of the interval. |
| `variant` | `'tiegcm'` (default, the reference) or `'paper'`. |
| `params` | dict that overrides pattern constants (see [Parameters](#the-parameters)). |

**Drivers.** The strength of the pattern is set by the cross-cap potential (CP). If it is not passed, it is derived from
Kp with TIE-GCM's empirical relation (next section), so the default is: **Kp from the DataManager -> CP -> pattern**.
By comes from the DataManager too.

**By is averaged by default.** By does not come from the store as the instantaneous value: the default is the **mean of
the previous 20 minutes of the 5-minute OMNI series** (`res='5min', avg=20`), the same default as `weimer05`.
`avg=None, res='1min'` gives the instantaneous 1-minute value; `avg`/`lag` give other windows. This is a choice
of this wrapper: TIE-GCM has no such rule, it uses whatever By the user supplies (a constant or a time series).
Values you pass in are used as given, with no averaging.

## In the survey

`survey.run_track` / `survey.run_grid` include the model (`models=["heelis82"]`, or the default of all applicable models). The
geographic samples are converted to AACGM with `aacgmv2`, as for `weimer05`, and the samples whose AACGM latitude is poleward of
30 deg are evaluated as aligned (time, mlat, mlt) samples; the others stay NaN (the model itself gives 0 there, but the survey marks
what a model does not cover as NaN), and so does everything above 2000 km (AACGM's limit). The variables are `heelis82_epot` [kV]
and the AACGM coordinates and drivers used, `heelis82_mlat`, `_mlt`, `_kp`, `_cp`, `_by`. `heelis_kw={...}` passes any
`models.heelis82` keyword:

```python
ds = survey.run_track(times, lats, lons, alts, models=["heelis82"])                     # Kp and By from the store
ds = survey.run_grid(times=[t], lats=np.arange(50, 90.1, 5.), lons=np.arange(-180, 180, 10.), alts=[400.],
                     models=["heelis82"], heelis_kw=dict(kp=5., by=-3., variant="paper"))
```

The geographic pre-filter is 15 deg (30 deg would miss the American sector, where AACGM 30 deg is near 20 deg geographic);
the drivers are looked up only if some sample is inside, so a low-latitude track needs nothing from the store.

## Kp -> cross-cap potential

The cross-cap potential (`ctpoten` in TIE-GCM, in kV) is the potential drop across the polar cap, between the morning and
the evening convection cells. When `ctpoten` is not given, TIE-GCM (sub `getgpi`, `gpi.F`) interpolates the 3-hourly
Kp to the model time and applies `ctpoten_from_kp` (`util.F`), which is ported unchanged, with its comments, to
`kp_ctpoten.f90`:

```
ctpoten = 15 + 15 Kp + 0.8 Kp^2      (kV; "formula given by Wenbin based on data fitting", LQIAN 2007; Kp in 0..9)
```

| Kp | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| CP [kV] | 15 | 30.8 | 48.2 | 67.2 | 87.8 | 110 | 133.8 | 159.2 | 186.2 | 214.8 |

The Kp interpolation is linear in time between the centres of the 3-h intervals (01:30, 04:30, ... UT), with
negative values set to zero, as in `gpi.F`. (A comment at the top of `gpi.F` still says `ctpoten = 29 + 11 Kp`; the code
uses the quadratic relation above, and so do we.) A Kp outside 0-9 stops a TIE-GCM run; here it gives NaN and a warning.

## What the function computes

The potential is the product of a latitude function and a local-time function, on a pattern of two convection cells
inside a circle of colatitude `theta0` (the *convection reversal boundary*) around a centre a little off the magnetic pole.

1. **Coordinates.** The point is rotated into circle coordinates: its colatitude `theta` from the pattern centre (offset
   `offc` toward midnight and `dskofc` toward dusk) and its local-time angle from noon. TIE-GCM's angle runs eastward
   (dusk at +90 deg); the paper's runs clockwise, so the two are mirror images that both give a positive morning cell.
2. **Local-time function F(phi, theta)** (eqs. 4-5 of the paper, `flwv32`). Around the dayside entrance `phid` and the
   nightside entrance `phin` the potential moves from the morning-cell value `psim` to the evening-cell value `psie`
   along a cosine, over a zone whose half-widths (`phidp0`, `phidm0`, `phinp0`, `phinm0`) grow with distance from the
   boundary as `phi + ((theta-theta0)/theta0)^2 (90 deg - phi)`, and stay `psim` / `psie` elsewhere. TIE-GCM uses 90 deg for
   all four half-widths, which makes the zones fill almost the whole circle (the "plateaus" of the paper's Fig. 4 mostly
   disappear).
3. **Latitude function.** In TIE-GCM (`variant='tiegcm'`):
   - outside the circle (`theta >= theta0`): `F * (sin(theta)/sin(theta0))^rr1 * exp(7 (1 - max(sin(theta), sin(theta_r1))/sin(theta_r1)))`,
     with `rr1 = -2.6` and `theta_r1 = theta0 + 11.3 deg` (the comment in the code: "average AMIE results show r1 = -2.6 for
     11.3 degrees beyond theta0");
   - inside: a cubic in `x = theta/theta0` that goes from the centre potential `pcen` (x = 0) to `F` at the boundary (x = 1) with a
     gradient set by `F` at the antipodal point, and zero slope at the boundary.
4. **Below 30 deg** the potential is set to 0 (TIE-GCM's `potm`); above it, values are in kV.

### The parameters

`theta` and `phi` as above. TIE-GCM computes all of the following in `aurora_cons` (`aurora.F`), once per time step,
from the two inputs `ctpoten` (CP, kV) and `byimf` (By, nT, limited to -11..+7); `heelis_cons` in `heelis82_modules.f90`
is that block. Values are northern hemisphere; the southern one has the sign of By reversed where By appears.

| Symbol | Meaning | Paper | TIE-GCM (variable) | `params` key |
|---|---|---|---|---|
| theta0 | colatitude of the convection reversal boundary: radius of the circle that holds the two cells, from the pattern centre | < 20 deg (15 deg in Fig. 4) | `theta0` = (-3.80 + 8.48 CP^0.1875) deg: 10.3 / 13.5 / 17.4 deg at 15 / 45 / 134 kV | `theta0` [deg] |
| psim | potential of the morning cell boundary (positive) | +20 kV (Fig. 4) | `psim` = +0.44 CP | `psim` [kV] |
| psie | potential of the evening cell boundary (negative) | -35 kV (Fig. 4) | `psie` = -0.56 CP (so psim - psie = CP) | `psie` [kV] |
| pcen | potential at the pattern centre | not in the paper (0 at the pole) | `pcen` = (-0.168 -/+ 0.027 By) CP | `pcen` [kV] |
| phid | centre of the dayside convergence zone ("throat"), where the potential switches from psim to psie | 15 deg clockwise from noon (Fig. 4) | `phid` = MLT 9.39 -/+ 0.21 By | `phid` [MLT h] |
| phin | centre of the nightside zone (Harang discontinuity) | 180 deg (Fig. 4) | `phin` = MLT 23.50 -/+ 0.15 By | `phin` [MLT h] |
| phi_d+, phi_d- | half-widths of the dayside zone at theta0 (morning / afternoon side) | 22.5 deg (Fig. 4); 20-60 deg (Fig. 5) | `phidp0`, `phidm0` = 90 deg | `phidp`, `phidm` [deg] |
| phi_n+, phi_n- | half-widths of the nightside zone (dusk / dawn side) | 55 deg (Fig. 4); 15-40 deg (Fig. 6) | `phinp0`, `phinm0` = 90 deg | `phinp`, `phinm` [deg] |
| offc, dskofc | offset of the pattern centre from the magnetic pole, toward midnight / toward dusk | 5 deg toward midnight (Figs. 5-6) | `offc` = 1.1 deg; `dskofc` = -0.08 -/+ 0.15 By deg | `offc`, `dskofc` [deg] |
| r1 | fall-off index equatorward of the boundary | -4 | `rr1` = -2.6 (and the exp cutoff) | `r1` |
| r2, theta_c | index and phase angle of the polar-cap function | 2 (or 1) and 14 deg | not used (cubic + `pcen`) | `r2`, `thetac` [deg] (paper) |
| theta1, theta2 | equatorward / poleward edges of the smoothing region around theta0 | free; 2 deg wide in Fig. 3 | not used | `dtheta1` = theta1 - theta0, `dtheta2` = theta0 - theta2 [deg] (paper) |

`params` overrides apply to both hemispheres. Only `ctpoten` and By drive TIE-GCM's values; everything else in the
table is a constant or a function of those two inside `aurora_cons`. In particular the paper has no dependence on
CP or By at all (its parameters are illustrative: "no attempt to match a data set"), so `theta0(CP)`, the 0.44 / -0.56
split, `pcen` and the By terms are TIE-GCM's, with comments crediting W. Wang (2008, the By terms and limits) and L. Qian
(2007, the Kp relation).

## Two variants

`variant='tiegcm'` (default) is the reference: the original expressions of `heelis.F`, identical to a TIE-GCM run.

`variant='paper'` replaces only the last step (the latitude function) by that of Heelis et al. (1982): `psi = G(theta) F(phi, theta)`.
G is the power law `A1 (sin theta/sin theta0)^r1` equatorward of `theta1`, elliptical arcs `sqrt(1 - (theta - theta0)^2/B)`
between `theta2` and `theta1`, and `A2 [(sin(theta + theta_c)/sin theta0)^r2 - (sin theta_c/sin theta0)^r2]` poleward of
`theta2`. The arcs remove the discontinuity of the north-south electric field at the boundary that a plain
`(sin theta/sin theta0)^r` has; TIE-GCM does not have them, and its field jumps by about 4 kV/deg across the boundary
(the paper calls that "an obvious discrepancy"). `A1, B1, A2, B2` make G and its slope continuous (closed form, in
`paper_cons`), `G(theta0) = 1` so the extremes are still `psim`/`psie`, and `G(0) = 0` (the pole has zero potential;
TIE-GCM's is about -9 kV at 45 kV). The defaults are the paper's "typical values", `r1 = -4`, `r2 = 2`, `theta_c = 14 deg`, and a
2-deg wide region, `theta1 = theta0 + 1`, `theta2 = theta0 - 1` (the paper does not say how its 2 deg is split).
Everything else (theta0, psim, psie, phid, phin, half-widths, offsets) is TIE-GCM's. **This is a hybrid**: the paper's function
driven by TIE-GCM's parameterization. Its constants were tuned together with `rr1 = -2.6`, the cutoff and `pcen`, and the
result has not been validated against data. Pass `params=` to use the paper's own numbers, for example its Fig. 4
(`psim=20, psie=-35, phid=11, phin=24, phidp=phidm=22.5, phinp=phinm=55, theta0=15, offc=0, dskofc=0, pcen=0`).

## Differences from a TIE-GCM run

The aim is the same potential; these are the places where it is not literally what TIE-GCM does:

- **It is the empirical potential, not TIE-GCM's final one.** TIE-GCM evaluates this pattern (`phihm`, sub `potm`) and
  combines it with its dynamo solution (`pdynamo.F`, using the fractional presence of the dynamo equation `pfrac` from
  `colath.F`; read from the source, not traced in a run). The function returns `phihm`; it does not contain the dynamo
  solution, nor the blending with it.
- **Coordinates.** TIE-GCM evaluates on its own 97 x 80 magnetic grid (a dipole/apex system), with the sun's magnetic
  longitude from its own routine. Here the caller gives magnetic latitude and MLT; with AACGM the result is close, not
  identical (the same approximation is made for `weimer05`).
- **Points at or below 30 deg** get 0, as in TIE-GCM. Points can be anywhere else, not only on TIE-GCM's grid.
- **By** defaults to a 20-minute mean of the OMNI 5-minute series (TIE-GCM uses the value the user supplies).
- **Kp out of range** gives NaN and a warning (TIE-GCM stops the run); the AMIE branch of `aurora_cons` is left out.
- Nothing here depends on the hemispheric power, which TIE-GCM only uses for the auroral precipitation.

## Files

| File | What it is |
|---|---|
| `heelis.F` | TIE-GCM 2.0's `heelis.F`, **modified** (header notice, dated): `flwv32` imports `hvariant` and branches to the paper's latitude function; `gpaper` is added. The rest is byte-for-byte the original. |
| `heelis82_modules.f90` | replacements for the four TIE-GCM modules `heelis.F` imports; `heelis_cons` is the part of `aurora_cons` (`aurora.F`) that sets the parameters, with its comments; `paper_cons` and the overrides are new |
| `kp_ctpoten.f90` | `ctpoten_from_kp` from `util.F`, comments included (returns NaN instead of stopping) |
| `heelis82_batch.f90` | the driver (ours): one call for many points, groups the points by drivers and hemisphere, calls `flwv32` |
| `heelis82f.pyf` | hand-written f2py signature; only the driver is exposed |
| `tiegcmlicense.txt` | the TIE-GCM license (see below); installed with the package |

The model needs 8-byte reals (`-fdefault-real-8`, as TIE-GCM is built with gfortran). There are no data files.

## Validation

`tests/test_heelis82.py`:

- `variant='tiegcm'` against the **untouched** `heelis.F`: the original file, with its parameter block transcribed
  independently from `aurora.F`, was run over 3600 points (3 CP x 6 By x both hemispheres x 10 latitudes x 10 MLT); the
  two agree to 2e-10 V (rounding). 41 of those values are in the tests.
- Kp -> CP against the formula, CP over Kp, the By limits, the north-south mirror (south = north with By reversed),
  the cross-cap potential (max - min of the potential = CP), the sign of the cells, cubes / aligned samples.
- `variant='paper'` against an independent implementation of G (closed form) times F (eqs. 4-5) with the paper's Fig. 4
  constants, at ten colatitudes from the pole to 9 deg beyond the boundary: within 0.04 V of a 55 kV swing (the test
  allows 0.5 V); and that G has no kink at the boundary where the TIE-GCM variant has one.
- The drivers: By averaged (20 min, 5-min series by default), Kp interpolated between interval centres like `gpi.F`.
- `examples/plot_heelis82_fig04.py` reproduces Fig. 4 of the paper (`variant='paper'` with the paper's constants: the
  +20 / -35 kV plateaus at the local times the paper gives, and the potential at the centres of the two convergence
  zones, checked to 0.1 kV), compares the two variants across the boundary, and maps the pattern for Kp = 3. It needs
  no network: Kp, By and the constants are given.

What has **not** been checked is the output of an actual TIE-GCM run (only its source), nor the paper's Figures 1-3
(the latitude profile of TIE-GCM is not the paper's, and the paper variant uses assumed smoothing widths).

## License

`heelis.F` and the parts derived from TIE-GCM are **NCAR TIE-GCM 2.0** software, under its *Open Source Academic Research
License Agreement* (`tiegcmlicense.txt`, which travels with them): use is for **research, academic and non-profit
purposes only**, the software **may not be used for operational purposes** (products for forecast, nowcast or hindcast of
the atmospheric state), it may not be sold, modifications must carry change notices with their dates (this folder's
files do), and derivative works are to stay open source. The reference the license asks for:

> This software is part of the NCAR TIE-GCM. Use is governed by the Open Source Academic Research License Agreement
> contained in the file tiegcmlicense.txt.

Section 2 of the license reserves the Weimer model in TIE-GCM (`wei05sc.F`) to its author; nothing of it is used here.
The academic "rules of the road" of NSF CEDAR and NASA TIMED apply to publications that use these results.
