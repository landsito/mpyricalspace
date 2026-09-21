# Weimer (2005) high-latitude electric potential and field-aligned current

Notes on running the Weimer model in mpyricalspace. The short version is in the
[main README](../../README.md#use);
this file has the options and limits, and, further down, how the code maps to the papers and how it compares with
the papers' figures.

> Weimer, D. R. (2005a), Improved ionospheric electrodynamic models and application to
> calculating Joule heating rates, J. Geophys. Res., 110, A05306,
> [doi:10.1029/2004JA010884](https://doi.org/10.1029/2004JA010884)
>
> Weimer, D. R. (2005b), Predicting surface geomagnetic variations using ionospheric
> electrodynamic models, J. Geophys. Res., 110, A12307,
> [doi:10.1029/2005JA011270](https://doi.org/10.1029/2005JA011270)

The bundled Fortran, translated by B. Foster (HAO/NCAR), is the spherical-cap-harmonic (SCHA) revision of the 2005a
models that 2005b describes (degrees up to 12, orders up to 2, cap-size-dependent degrees, one shared boundary). It
gives the **electric potential** and the **field-aligned current** poleward of the model's low-latitude boundary,
driven by the IMF, the solar wind and the dipole tilt. `Empirical.run_weimer05` is a thin wrapper over
`models.weimer05` (same arguments, `obj.time` as the times).

**Contents:** [Inputs, coordinates and output](#inputs-coordinates-and-output) ·
[Driver averaging](#driver-averaging) · [Examples](#examples) · [In the survey](#in-the-survey) ·
[How the code maps to the papers](#what-the-papers-say-and-where-it-lives-in-the-code) ·
[Validation](#validation-against-the-papers-figures) · [Build notes](#build-notes) · [Credits and license](#credits-and-license)

## Inputs, coordinates and output

- **Inputs.** `time`; `mlat`, `mlt` (below); `by`, `bz` [nT, GSM], `vsw` [km/s], `nsw` [cm⁻³] as a scalar or one
  value per time — omitted, they come from `DataManager`; `tilt` [deg] — omitted, `models.dipole_tilt(time)` computes
  it (IGRF-14 dipole, checked against astropy's Sun to 0.005°); `res`, `avg`, `lag` — how store-provided drivers are
  averaged (defaults `res='5min', avg=20, lag=0`; see [Driver averaging](#driver-averaging)).
- **Coordinates.** The model was built in **AACGM** (altitude-adjusted corrected geomagnetic) latitude and magnetic
  local time (2005a, section 2 and appendix A). When you call it directly, converting from geographic coordinates is up
  to you (`aacgmv2`, installed with the package, does it; the survey uses it automatically). The paper notes that
  modified-apex and AACGM coordinates differ very little at high latitudes, so apex-type latitudes such as the
  quasi-dipole `mlat` from `run_hltwim` should be close there; that has not been quantified here. `fac` is a density on
  the magnetic grid: mapping it onto a geographic grid needs a compensation for the non-uniform mapping between the two
  systems (Richmond, 1995, as the papers note).
- **Output.** `epot` in **kV** and `fac` in **µA/m² at 110 km, positive downward** (into the ionosphere), as in the
  paper's figures, plus the drivers used along `time`. NaN equatorward of the cap, and wherever a driver is
  missing. **Southern hemisphere:** pass negative latitudes; they are handled as the paper prescribes (2005a,
  section 5) — the northern model mirrored, with the signs of By and of the tilt reversed.
  A cap size outside what the tables support (roughly 8°–66° of colatitude; only extreme storms leave it) gives
  NaN and a `RuntimeWarning` instead of crashing. Normal conditions give a boundary of about 24°–36° colatitude.
- **Range of validity.** The fits used IMF magnitudes below roughly 15 nT; beyond that the nonlinear saturation
  law extrapolates (the paper's example: 250 kV at Bz = −50 nT) and it advises caution.
- **Speed.** All points of a time step share one model set-up: a full day of 1-min steps on a 97 × 25 grid
  takes about 1.5 s.

## Driver averaging

Drivers you do not pass (`by`, `bz`, `vsw`, `nsw`) come from the index store (NOAA OMNI; see
[Geophysical indices](../../README.md#geophysical-indices)). **By default each is the mean of the previous 20 minutes
of the 5-minute series** (`res='5min', avg=20`), which is how Weimer (2005b) drives the model — not the instantaneous
value.

| Call | Drivers taken from the store |
|---|---|
| default (`res='5min', avg=20`) | mean of the previous 20 min of the 5-min series (Weimer 2005b) |
| `res='1min', avg=None` | the instantaneous 1-min value |
| `avg=45, lag=10` | the recipe the coefficients were fitted with: the propagated IMF from 55 to 10 min before each step (Weimer 2005a) |
| values you pass yourself | used as given, never averaged |

OMNI is already time-shifted to the bow shock, whereas the papers propagate to the magnetopause nose (the difference is
small against the averaging window). OMNI speed and density have gaps (10–13 May 2024, for example, has all four
drivers on about two thirds of the minutes), which give NaN rows unless you pass `nsw=` yourself; the papers
substitute average values (2005a) or hold the last known density (2005b, which notes it matters little).

## Examples

```python
import numpy as np
from datetime import datetime, timedelta
from mpyricalspace import Predictor, models

# one time step on a magnetic grid: a (mlat, mlt) cube, drivers given by hand
ds = models.weimer05(datetime(2024, 5, 11, 12), mlat=np.arange(50., 90.1, 1.), mlt=np.arange(0., 24., 1.),
                     by=0., bz=-5., vsw=450., nsw=9.)
ds.epot            # kV      (mlat, mlt)
ds.fac             # uA/m^2  (mlat, mlt), positive = downward, at 110 km

# through the facade: drivers from the index store (default: the mean of the previous 20 min of the 5-min series)
obj = Predictor.Empirical(); obj.set_period(d0, dn, freq=timedelta(minutes=5))
ds = obj.run_weimer05(mlat, mlt)                                # (time, mlat, mlt)
ds = obj.run_weimer05(mlat, mlt, res='1min', avg=45, lag=10)    # the recipe the coefficients were fitted with

# a satellite track: equal-length time / mlat / mlt arrays are aligned samples
ds = models.weimer05(times, track_mlat, track_mlt)
```

`examples/plot_weimer05_fig02.py` reproduces the paper's Figure 2 (nine polar maps of the potential) and checks the
extremes against the values printed in it; it runs as part of `tests/test_examples.py`.

## In the survey

`survey.run_track` / `survey.run_grid` include the model (`models=["weimer05"]`; the default `models=None` run
includes it too). The geographic samples are converted to AACGM with
[`aacgmv2`](https://pypi.org/project/aacgmv2/) (installed with the package). Only samples poleward of 30° geographic
latitude are converted, and nothing above 2000 km (AACGM's limit) is evaluated; everything else is `NaN`. It returns
`weimer05_epot` [kV] and `weimer05_fac` [µA/m², + downward] plus the AACGM coordinates and the drivers it used
(`weimer05_mlat`, `_mlt`, `_by`, `_bz`, `_vsw`, `_nsw`, `_tilt`). If `aacgmv2` is missing anyway (for example after a
failed build), the model lands in `ds.attrs["skipped"]` with an explanatory message, and only if some sample actually
needs it. Drivers come from the index store as above; `weimer_kw={...}` passes any `models.weimer05` keyword.

```python
ds = survey.run_track(times, lats, lons, alts, models=["weimer05"])                     # drivers from the store
ds = survey.run_track(times, lats, lons, alts, models=["weimer05"],
                      weimer_kw=dict(by=0., bz=-5., vsw=450., nsw=9., tilt=0.))         # or given by hand
ds = survey.run_grid(times=[t], lats=np.arange(60, 90.1, 5.), lons=np.arange(-180, 180, 10.),
                     alts=[400.], models=["weimer05"])                                  # a geographic grid
```

`fac` is the current density at 110 km, not at the sample's altitude, so along a satellite track only `epot`
(constant along field lines) is directly comparable to what the satellite sees. For maps on a *magnetic* grid use
`run_weimer05` directly.

## What the papers say, and where it lives in the code

The equations of 2005a (boundary regression, eq. 1–3; coefficient response, eq. 6–7) were checked line by line
against `setmodel` / `setboundary` and the data files, and match:

| Paper | In mpyricalspace |
|---|---|
| boundary radius R = Σ Bᵢ Wᵢ, W = [1, cos θ, E(B_T), E(B_T) cos θ, V_SW, P_SW, (AL)], a circle offset **4.2° toward 0 MLT**; p₁, p₂ = −0.16, 0.25 (2005a eq. 1–3) | `setboundary` → `bndyfitr`; `xc = 4.2`, `x = [1, cosθ, btx, btx·cosθ, v, P]` (no AL); `W05scBndy.dat` holds the two coefficient sets (with and without AL) and 0.1595, 0.2478 |
| Aₘ = Σₙ Σᵢ Cᵢₙₘ Xᵢ fₙ(θ), X = [1, E(B_T V_SW) with B_T V_SW as the solar-wind electric field in mV/m, sin t, sin² t, P_SW, (AL)]; **15 coefficients per basis function without AL**; p₁, p₂ = −1.33, 0.47 (electric), −1.76, 0.70 (magnetic) (eq. 6–7, sec. 4) | `setmodel`: `swe = (1−exp(−sw·p₁))·sw^p₂` with `sw = Bt·V/1000`, `swp = n·V²·1.6726e-6`, the 15-term vector `a(1:15)`, `d1_pot = 15`; `ex_pot` = (0.47, 1.33) / (0.70, 1.76) in the data files |
| θ = the IMF clock angle in the **GSM** Y–Z plane, 0° = northward | `angle = atan2(By, Bz)` |
| potential from SCHA functions Pₙₖ₍ₘ₎ᵐ(cos θ)(g cos mφ + h sin mφ), k ≤ 12, m ≤ 2, degree *n* depending on the cap half-angle θ₀; k−m odd terms only, so the potential is **zero on the boundary** (2005b eq. 11) | `read_data` (`ab`, `ls`, `ms`); `nkmlookup` interpolates *n* in `SCHAtable.dat`; `scplm` caches the basis functions per θ₀ |
| electric and magnetic potentials share one boundary | one `bndyfitr` for both |
| J∥ = ∇²ₛψ / μ₀ (eq. 12); the figures plot **downward positive**, at 110 km | `mpfac`: the same series weighted by n(n+1) (the surface-Laplacian eigenvalue), scaled to the 110 km shell — hence `setmodel('bpot')` before `mpfac` |
| southern hemisphere: reverse the signs of By and of the tilt (2005a sec. 5) | done in `weimer05()` |
| models were fitted with the propagated IMF averaged from 55 to 10 min before each pass (45 min; 2005a sec. 2); 2005b drives them every 5 min with the previous 20 min | `avg=45, lag=10` / `res='5min', avg=20` on OMNI (already shifted to the bow shock; the paper says correcting for the bow-shock position is negligible against a 45-min average) |
| missing plasma data: average values substituted (2005a); the last known density held (2005b) | not automatic — pass `nsw=` |

Two implementation details that are *not* in the papers: below 1 nT of IMF the code removes the clock-angle
dependence (`bt < 1` branches in `setmodel` and `setboundary`), and `weimer05_batch` returns NaN outside the
tables' cap-size range.

**Not part of this Fortran, though the papers use them:** the optional AL-index (substorm) component of the
models; the electric field (−∇ of the potential), Joule heating / Poynting flux (2005a sec. 6), the equivalent
current and the ground-level ΔB (2005b: Biot–Savart at 110 km, Hall/Pedersen ratio 1.5).
Caveats stated in the papers: the FAC has **artifacts near the geomagnetic pole** (second derivatives of the
highest-degree terms); the models do not capture the response to sudden IMF or dynamic-pressure changes, nor the
sub-auroral penetration field; and far-out-of-range IMF is extrapolation (above).

## Validation against the paper's figures

2005a Figures 2 and 4 give the potential and FAC extremes for BT = 5 nT at nine clock angles (V = 450 km/s,
N = 4 cm⁻³, tilt 0), and the total downward current for each. Those figures show the earlier *hybrid* version, so
they are a benchmark, not an exact target. `examples/plot_weimer05_fig02.py` reproduces Figure 2 and fails if
any extreme drifts by more than 4 kV from the printed values:

| IMF clock angle | electric potential min / max [kV], paper → here | total downward current [MA], paper → here |
|---|---|---|
| 0° (north) | −12 / 12 → −13.1 / 13.0 | 1.4 → 1.16 |
| 90° | −40 / 29 → −37.1 / 30.2 | 2.3 → 2.26 |
| 180° (south) | −55 / 47 → −54.7 / 47.5 | 4.1 → 3.73 |
| 270° | −30 / 28 → −30.9 / 29.2 | 2.2 → 2.21 |

Potentials agree to about 4% (median over all nine cases). The **FAC peak values are lower** than the paper's
(for 180°: −0.80 / +0.69 µA/m² here against ±1.03) although the *totals* and the sign pattern (region 1
downward at dawn, upward at dusk) agree; 2005a itself notes that the hybrid model's FAC densities are
"noticeably higher because of the larger second derivatives", and 2005b that the SCHA maps are smoother.

## Build notes

The model is built by meson together with the rest of the package (`src/meson.build`). Its compiler flags are
requirements of the model, not tuning: `-fdefault-real-8` (default `real` and unsuffixed literals need 8 bytes) and
`-fno-automatic` (locals such as `skip` are read before they are set and rely on zero-initialised static storage).
Only the driver `weimer05_batch` is exposed to python; the model's own routines are compiled in and called from it. The `.dat` files are read from `_data/weimer05/` (an installed
package) or `src/weimer05/` (editable install).

## Credits and license

The Fortran in this folder is a Fortran 90 translation, by B. Foster (HAO/NCAR, 2008), of D. Weimer's IDL model. It
ships with its own [`LICENSE`](LICENSE), **separate from this package's MIT license**: non-commercial scientific
research and education, with attribution to the author and publications retained; commercial use or redistribution
needs separate permission. The file is installed next to the data (`_data/weimer05/LICENSE`) and must stay with any
redistribution. Please cite Weimer (2005a, 2005b) and credit the translation.
