# mPyricalspace

mPyricalspace provides empirical predictions of the near-space environment: it
defines the geophysical state of the ionosphere and thermosphere from the
climatological empirical models available.

## Installation

A pip-installable release is planned; for now install from a source checkout.
Build-time system requirements (the models are Fortran/C extensions):

- a Fortran compiler — `gfortran` (`brew install gcc` / `apt-get install gfortran`)
- **GSL** — used by the JVDM / EEJ / EEF (Alken) models:
  `brew install gsl` / `apt-get install libgsl-dev` / `conda install -c conda-forge gsl`

Then:

```bash
git clone https://github.com/landsito/mpyricalspace.git
cd mpyricalspace
pip install .
```

The Python dependencies (numpy, pandas, xarray, pymsis, duckdb, pymongo, pillow,
requests) install automatically. Model coefficient files (HWM, IRI, EEJ/EEF, ACE)
and P. Alken's `ndlinear` are bundled — nothing else to *need* downloading. HWM14/
HWM07, IRI, IGRF and EEJM1/EEJM2 are the exception: at build time each tries to
verify its own vendored source against map.nrl.navy.mil / irimodel.org /
ngdc.noaa.gov / geomag.colorado.edu and refresh anything that changed (see
`src/hwm/README.md` / `src/iri/README.md` / `src/igrf/README.md` /
`src/eejm/README.md`); this needs network access but never blocks the build — no
connectivity just means the bundled copy is used as-is. The index backend
(F10.7/Kp/Ap/AE, separate from this) is chosen at *runtime*, not at install time.

Check the build with `python -m mpyricalspace doctor` — a table with one row
per compiled Fortran extension (a partial build leaves some missing). See
[Command line](#command-line) below for the full set of subcommands.

### Editable install (for development)

Same as above, but with `-e` so the package tracks your working tree (the
extensions are rebuilt automatically on import when a source file changes):

```bash
git clone https://github.com/landsito/mpyricalspace.git
cd mpyricalspace
pip install -e . --no-build-isolation
```

### macOS + conda notes

- If `pkg-config` can't find GSL (conda's `pkg-config` doesn't look in Homebrew's
  path), point it there:

  ```bash
  export PKG_CONFIG_PATH=/opt/homebrew/lib/pkgconfig      # Apple Silicon
  ```

- If linking the static Fortran libs fails with `ar` / `ranlib` errors (conda's
  binutils vs Apple's linker), force Apple's toolchain with the bundled native file:

  ```bash
  pip install . --config-settings=setup-args=--native-file=$PWD/native.txt
  ```

  or just `export AR=/usr/bin/ar RANLIB=/usr/bin/ranlib` before `pip install`.

- If linking fails with `ld: warning: ignoring file .../libSystem.tbd, malformed file`
  (`error: unknown architecture arm64e.x1-macos`) followed by dozens of
  `Undefined symbols` (`_sin`, `_malloc`, ...), conda's `ld64` is older than the newest
  Command Line Tools SDK (`MacOSX.sdk` may point at a release it cannot parse). Nothing is
  wrong with the code or the architecture (both are arm64) — point the build at an older
  SDK that is installed, once per environment:

  ```bash
  ls /Library/Developer/CommandLineTools/SDKs/
  conda env config vars set SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.5.sdk -n <env>
  conda deactivate && conda activate <env>       # the variable is set on activation
  ```

## Use

```python
from datetime import datetime, timedelta
from mpyricalspace.Predictor import Empirical

obj = Empirical()
obj.set_period(datetime(2024, 5, 11), datetime(2024, 5, 12), freq=timedelta(minutes=15))
obj.set_location(-11.958256, -76.859012, 250)      # hwm / msis / igrf need this
ds = obj.run_rocsat()                              # every run_* returns a fresh xr.Dataset
ds = obj.run_scherliessfejer()                     # quiet + storm vertical drift
ds = obj.run_jvdm1()                               # Alken JULIA vertical drift, 150 km (08-16 LT)
ds = obj.run_hwm(version=2014)                     # or version=2007 / 1993
ds = obj.run_hltwim(lat=70., lon=-100.)            # high-latitude winds (|MLAT|>40)
ds = obj.run_eej(version=2)                        # equatorial electrojet; or version=1
ds = obj.run_eej(version=2, model="oersted")       # ... from CHAMP (default) / Orsted / SAC-C
ds = obj.run_eef()                                 # equatorial electric field
ds = obj.run_msis(alt=300.)                        # NRLMSIS 2.1 thermosphere (via pymsis)
ds = obj.run_iri(alts=300., version=2026)          # or version=2020 / 2016 / 2012 / 2007 / 2001
ds = obj.run_igrf(version=14)                      # or version=13 / 12 / 11 / 10 / 9
ds = obj.run_ppeefm1()                             # Manoj & Maus PPEF web service (any date)
ds = obj.run_ppeefm1(nativelypackaged_code=True)   # same model, run locally from bundled ACE (2001-2007 only)
```

`Empirical` is a thin facade; the models are pure functions in
`mpyricalspace.models` (`rocsat_drift`, `scherliess_fejer`,
`jvdm1_drift`, `eej`, `eef`, `hwm`, `hltwim`, `msis`, `igrf`,
`iri`) and can be called directly. The shared n-d grid builder / cwd guard /
Dataset packer live in `mpyricalspace._grid`.

### Many models at once — along a track or over a grid

`mpyricalspace.survey` runs every *applicable* model for a set of points and
merges the results into one `xr.Dataset`; each model is evaluated only where it
is defined, everything else is `NaN`.

```python
from mpyricalspace import survey

# a satellite ephemeris: equal-length time / lat / lon / alt(km) arrays
ds = survey.run_track(times, lats, lons, alts)     # models=None -> every applicable model
ds.attrs["models"]                                 # what ran; ds.attrs["skipped"] what didn't

ds = survey.run_track(times, lats, lons, alts, models=["iri", "msis"])   # or a subset

# or a meshgrid
ds = survey.run_grid(times=[t], lats=np.arange(-60, 61, 5.),
                     lons=np.arange(-180, 180, 5.), alts=[350.])

list(survey.model_domains())                       # the allowed model names
survey.model_domains()                             # + kind / altitude window / mag-lat cap each
```

The high-latitude models (`hltwim`, `weimer05`) are evaluated only poleward of their boundary. **`weimer05`**
converts the geographic samples to AACGM with [`aacgmv2`](https://pypi.org/project/aacgmv2/)
(installed with the package), evaluates them, and returns `weimer05_epot` [kV], `weimer05_fac`
[µA/m², + downward] plus the AACGM coordinates and the drivers it used (`weimer05_mlat`, `_mlt`, `_by`, `_bz`,
`_vsw`, `_nsw`, `_tilt`). Only samples poleward of 30° geographic latitude are converted, and nothing above 2000 km
(AACGM's limit) is evaluated; everything else is `NaN`. If `aacgmv2` is missing anyway (for example after a failed build), the model lands in `ds.attrs["skipped"]` with
an explanatory message, and only if some sample actually needs it. By default the IMF / solar wind come from the
index store as the mean of the previous 20 min of the 5-min OMNI series (see
[`run_weimer05`](#empiricalrun_weimer05--weimer-high-latitude-electric-potential-and-field-aligned-current));
`weimer_kw={...}` passes any `models.weimer05` keyword, e.g. explicit drivers:

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

The global models (`igrf`, `msis`, `hwm`, `iri`) run at every sample. The
equatorial-electrodynamics models describe a magnetic-equator parameter vs
longitude and local time, so each is masked to the dip equator (IGRF dip
latitude, capped per model) *and* to its altitude regime:

- `sf` — Scherliess-Fejer *quiet* vertical drift; F region 200-900 km; dip lat within 2.5 deg
- `rocsat` — ROCSAT-1 quiet vertical drift; F region; within 3 deg
- `eej` — Alken equatorial electrojet current; E region 90-130 km; within 5 deg
- `eef` — Alken equatorial zonal electric field; 90-1000 km; within 10 deg
- `manoj` — Manoj-Maus RTEEF prompt-penetration equatorial E-field [mV/m]; F region 200-900 km; within 5 deg

Variables carry the model name (`igrf_B`, `hwm_u`, `iri_NMF2`, `eej`,
`sf_qvdrift`, `manoj_ppef` ...). `obj.run_track(lats, lons, alts)` /
`obj.run_grid(...)` are the same thing keyed off `obj.time`.

`manoj` is the local RTEEF, so it is ACE-limited to one calendar year in
2001-2007 (skipped, with a reason, otherwise). `sf` and `manoj` are two
*different* prompt-penetration models (storm climatology giving a drift [m/s] vs a
solar-wind transfer function giving the equatorial E-field [mV/m]) — pick one.
The Scherliess-Fejer *storm*
decomposition and `jvdm1` are Jicamarca / Peruvian-sector
models — use their own `run_*`.

The `examples/` directory has one runnable script per model (most reproduce a
figure from the model's own paper), plus `plot_survey_track.py` for the
multi-model `survey` runner. Each script's docstring explains what it does and
what it reproduces — read the header, then run it. Running them all is also a
quick check that a fresh build works:

```bash
python examples/plot_hwm14_fig03.py
pytest -q                            # or just the headless smoke test: pytest -q tests/test_examples.py
```

## Command line

Every subcommand needs a working install — i.e. `import mpyricalspace` must
succeed. If it does not, start with `doctor`.

```bash
python -m mpyricalspace doctor                 # table: does each compiled Fortran extension load?
python -m mpyricalspace doctor --json          # the same, machine-readable

python -m mpyricalspace status                 # local index-store coverage
python -m mpyricalspace fetch --from 2024-01-01 --to 2024-06-01 [--kpap-source ...] [--ae-source ...] [--include kpap ae dst asy sw1 sw5]
python -m mpyricalspace update [--since 2024-01-01]   # extend every series toward today
python -m mpyricalspace config --mongo-uri ... --mongo-user ... --data-dir ...

python -m mpyricalspace vscode                 # (re)install the companion VS Code extension
python -m mpyricalspace vscode status          # is it installed / was it auto-tried?
python -m mpyricalspace vscode build DIR       # just write the .vsix into DIR
```

`doctor` is also exposed as `mpyricalspace.build_info()` -> `{extension:
None if OK else "<error>"}` and `mpyricalspace.build_report()` (the table as a
string). `status` / `fetch` / `update` / `config` are covered under
[Geophysical indices](#geophysical-indices); `vscode` under [Companion VS Code extension](#companion-vs-code-extension).

## Companion VS Code extension

`mpyricalspace/vscode_ext/` is a small, **optional** VS Code extension. It runs
`doctor --json` and tints each `src/<model>/` folder in the Explorer green
(the model's extension compiled and loads) or red (missing / import error, with
the error in the tooltip); a status-bar item shows the `n/m` count. It is
editor-only — it writes nothing and touches neither the package nor the build.

Because Python wheels have no post-install hook, installation is *attempted*,
once, best-effort: on the first `python -m mpyricalspace ...` call and on the
first `import mpyricalspace` **from a VS Code terminal**. Failures are appended
to `~/.cache/mpyricalspace/vscode-ext.log` and never surface; a sentinel stops
it retrying. Install (or reinstall) it by hand with `python -m mpyricalspace
vscode`. Fully quit and reopen VS Code the first time. Opt out of the
auto-attempt with `MPYRICALSPACE_NO_VSCODE=1`; uninstall with
`code --uninstall-extension mpyricalspace.mpyricalspace-build-status`.

## Geophysical indices

Most models need F10.7 / Kp / Ap / AE. `mpyricalspace.DataManager` supplies them
transparently — you normally never touch it.

**Default — a local DuckDB + Parquet store.** It lives under
`$MPYRICALSPACE_DATA_DIR` (default `~/.cache/mpyricalspace`) and **auto-fills**
whatever date span a model asks for, downloading from the public sources below.
Drive it explicitly if you want to pre-fetch or inspect it:

```bash
python -m mpyricalspace status
python -m mpyricalspace fetch  --from 2024-01-01 --to 2024-06-01 [--kpap-source ...] [--ae-source ...] [--include kpap ae dst asy sw1 sw5]
python -m mpyricalspace update                       # extend every series toward today
```

**Optional — your own MongoDB.** If you already run a MongoDB with the
per-minute index schema (collection `geomagnetic_indices.resolution_1min`),
point `DataManager` at it and it is used instead of the local store. Give the
URL and credentials once and they are saved to `~/.config/mpyricalspace/config.json`:

```python
from mpyricalspace.DataManager import configure
configure(mongo_uri="mongodb://host:27017/", mongo_user="u", mongo_password="p")
```

```bash
python -m mpyricalspace config --mongo-uri mongodb://host:27017/ --mongo-user u --mongo-password p
```

or per-session via `$MPYRICALSPACE_MONGO_URI` (with `$MPYRICALSPACE_MONGO_USER`
/ `$MPYRICALSPACE_MONGO_PASSWORD`). `local=True` or `$MPYRICALSPACE_LOCAL=1`
ignores any configured MongoDB and forces the local store.

The same `fetch` / `update` can also **populate** a MongoDB (forward-filled
1-min docs), but only when writes are explicitly enabled, so a read-only database
is never touched by accident:

```bash
export MPYRICALSPACE_MONGO_WRITE=1        # or: config --mongo-writes, or DataManager(writes=True)
python -m mpyricalspace update
```

Download sources:

- **kp / ap / f107 / f107a** — GFZ Potsdam
  (`--kpap-source` potsdam | potsdam-legacy | noaa | canada | lisird)
- **ae / ao / al / au** — Kyoto WDC aeasy-cgi -> realtime-image digitisation -> NOAA OMNI
  (`--ae-source` auto | kyoto | images | omni)
- **dst** — Kyoto WDC
- **symh / asy** — Kyoto WDC aeasy-cgi
- **bt bx by bz vsw nsw tsw psw esw** — NOAA OMNIWeb 1-min IMF and solar wind
  (|B|, Bx, By, Bz in GSM [nT], flow speed [km/s], proton density [cm⁻³] and temperature [K],
  flow pressure [nPa], electric field [mV/m]); append `_5min` for the 5-min series
  (`by_5min`, `vsw_5min`, ...). 1-min starts in 1995, 5-min in 1981, and OMNIWeb lags real time by a few weeks,
  so the newest days come back NaN until it catches up. Groups `sw1` / `sw5` for `--include`.
  On the mongo backend these live in database `solar_wind`, collections `resolution_1min` / `resolution_5min`,
  under the long OMNI field names (e.g. `BY,_nT_(GSM)_(NOAA-OMNI)`); the short names above are translated on read
  and write, so a database that already follows that schema is used as is.
  ```python
  dm.get(dates, ['by', 'bz', 'vsw', 'nsw'])                         # 1-min
  dm.get_history(dates, 'bz', [-19, 0], timedelta(minutes=1))       # trailing 20 min of Bz
  ```

**Driving the Weimer model from these series.** The Weimer (2005a) coefficients were fitted with the propagated
IMF **averaged from 55 to 10 min before** each observation, and Weimer (2005b) drives the model every 5 min with
the mean of the **previous 20 min**. `run_weimer05(..., res='1min', avg=45, lag=10)` and `run_weimer05(..., res='5min', avg=20)` do those from this store —
the second is the **default** (see
[Weimer high-latitude electric potential](#empiricalrun_weimer05--weimer-high-latitude-electric-potential-and-field-aligned-current)).
Two things to keep in mind: OMNI is already time-shifted to the bow shock, whereas the papers propagate to the
magnetopause nose (the difference is small against the averaging window); and OMNI speed / density have gaps
(10–13 May 2024, for example, has all four drivers on about two thirds of the minutes), which give NaN rows unless
you pass `nsw=` yourself — the papers substitute averages or hold the last known density.

`run_scherliessfejer` evaluates the quiet + storm drift at every `obj.time`
and needs 28 h of AE before the earliest one (15-min AE, hourly as a fall-back);
gaps are filled with the last-known value. Pass `indices=` — an
`(n, >=2)` array of `[datetime, AE(, F10.7)]` — to bypass the store.

**`nativelypackaged_indices`**: a few models vendor their own index file as an
alternative to DataManager (most don't — they only ever use DataManager, there is
nothing else to switch). Every one of those follows the same rule: `False` (the
default, everywhere) always means DataManager; `True` means the model's own bundled
file instead. An explicit index passed by hand (`F107=`, `ap=`, `indices=` ...)
overrides either way.

- **IRI** (`run_iri`) ships `ig_rz.dat` (Rz12 / IG12) and `apf107.dat` (F10.7).
  `nativelypackaged_indices` only switches F10.7 / F10.7a (DataManager vs.
  `apf107.dat`, given explicitly with `F107=` / `F107a=` either way). `Rz12` /
  `IG12` have no DataManager equivalent and are **always** read from `ig_rz.dat` regardless of this flag; that file
  is copied on first use to `$MPYRICALSPACE_DATA_DIR/iri{12,16,20,26,07,01}`. Whenever
  a requested date runs past the maximum date in the file-cache, it gets refreshed from the CHAIN / ECHAIM mirror. This refresh never goes into DataManager's own store as it is never fed
  from ECHAIM. `python -m mpyricalspace update` refreshes
  `ig_rz.dat` / `apf107.dat` for every IRI version, or call
  `mpyricalspace.models.refresh_iri_indices()`.
- **MSIS** (`run_msis`) switches F10.7 / F10.7a / ap between DataManager and
  `pymsis`'s own cached index file (`pymsis.utils.get_f107_ap`, source: CelesTrak
  — entirely outside DataManager / ECHAIM).
- **Manoj & Maus local RTEEF** (`run_ppeefm1(nativelypackaged_code=True)` or `ief=`) switches
  the F10.7 fed to its climatological term between DataManager and the bundled SPIDR `F107.txt`; the
  prompt-penetration term never reads F10.7, so this flag leaves it unchanged.

`run_eej` / `run_eef` take longitude, local time and season from `time` /
`lon`, and F10.7 (used as the F10.7 / EUVAC proxy) from the store. The lunar
tide's dependence is built into the models (EEJ v2 and EEF); only the *lunar
local time* fed to them is approximated — a first-order formula that advances the
Moon-Sun elongation at the constant mean synodic rate, so the value can be
~20-25 min off. That is negligible against the ~12.4 h tidal period for
climatology; for lunar-tide studies pass `lunars=` computed from a real lunar
ephemeris. Points outside the model's local-time window (5-19 h; 7-17 h for EEF)
come back as `NaN`.

`run_eej(version=2, model=...)` picks which satellite the coefficients were fit
from: `"champ"` (default), `"oersted"` or `"sac-c"` (Alken & Maus 2007).
SAC-C flew a fixed local time, so its J has no local-time dependence.

`run_ppeefm1` (`nativelypackaged_code=False`, the default) calls the Manoj & Maus
PPEF web service (pass `data=` to use your own cached array).
`run_ppeefm1(nativelypackaged_code=True)` instead runs the *same* model (RTEEF —
Manoj, Maus & Alken) locally, from the ACE
solar-wind data bundled with the package: no network, fully reproducible, but valid
**only for 2001-2007**, within one calendar year, for spans of at most 5 days. The
ACE files ship gzip'd and are unpacked into `$MPYRICALSPACE_DATA_DIR` on first
use. Its climatological term is zero outside 07-17 h local time, so the drift
scaling is only meaningful during the day; expect close but not identical agreement
with the web service (which uses live, occasionally revised, ACE data).

`run_ppeefm1(ief=<array>)` feeds a **synthetic IEF Ey** [mV/m] straight into
the model's prompt-penetration filter (`TF.COF`), bypassing the ACE buffer —
useful for characterising the transfer function (see `examples/plot_mm_fig02.py`).
The array is 5-min cadence, length `(dn - d0)/300 s + 13` (one hour of priming
before `obj.time[0]`). The return is the raw model output in mV/m (`qef`,
`ppef`, `tef`) rather than the drift-scaled fields; `qef` is 0. See
`nativelypackaged_indices` above for where its F10.7 baseline comes from.

## Empirical models supported

Freely available to the community within freely available larger models 
or from development sites, or developed and provided by different researchers.

#### `Empirical.run_hwm` (`version=2014` (default) / `2007` / `1993`)

> Drob, D. P., et al. (2015), An update to the Horizontal Wind Model (HWM): The
> quiet time thermosphere, Earth and Space Science, 2, 301-319,
> [doi:10.1002/2014EA000089](https://doi.org/10.1002/2014EA000089)
>
> Emmert, J. T., et al. (2008), DWM07 global empirical model of upper thermospheric
> storm-induced disturbance winds, J. Geophys. Res., 113, A11319,
> [doi:10.1029/2008JA013541](https://doi.org/10.1029/2008JA013541)

#### `Empirical.run_hltwim`

> Dhadly, M. S., Emmert, J. T., Drob, D. P., McCormack, J. P., and Niciejewski,
> R. J. (2019), HL-TWiM empirical model of high-latitude upper thermospheric
> winds, J. Geophys. Res. Space Physics, 124, 6822-6841,
> [doi:10.1029/2019JA026779](https://doi.org/10.1029/2019JA026779)

#### `Empirical.run_iri` (`version=2026` (default) / `2020` / `2016` / `2012` / `2007` / `2001`)

> Bilitza, D., Truhlik, V., Yoshihara, O., and Moldwin, M. B. (2024),
> Development and improvement of the International Reference Ionosphere with
> special emphasis on the topside and extension to the plasmasphere,
> Ann. Geophys., 67(4), SA443,
> [doi:10.4401/ag-9145](https://doi.org/10.4401/ag-9145)

#### `Empirical.run_igrf` (`version='latest'` (default: the newest generation bundled, currently 14) / `14` / `13` / `12` / `11` / `10` / `9`)

> Beggan, C. D., Kloss, C., Amblard, P. et al. (2026),
> International geomagnetic reference field: the fourteenth generation,
> Earth Planets Space, 78, 127,
> [doi:10.1186/s40623-025-02360-0](https://doi.org/10.1186/s40623-025-02360-0)

#### `Empirical.run_jvdm1`

> Alken, P. (2009), A quiet time empirical model of equatorial vertical plasma
> drift in the Peruvian sector based on 150 km echoes,
> J. Geophys. Res., 114, A02308,
> [doi:10.1029/2008JA013751](https://doi.org/10.1029/2008JA013751)

#### `Empirical.run_eej` / `run_eef` (P. Alken's climatological EEJ and EEF models)

> EEJ — Alken, P., and S. Maus (2007), Spatio-temporal characterization of the
> equatorial electrojet from CHAMP, Orsted, and SAC-C satellite magnetic
> measurements, J. Geophys. Res., 112, A09305,
> [doi:10.1029/2007JA012524](https://doi.org/10.1029/2007JA012524)
>
> EEF — Alken, P., and S. Maus (2010), Electric fields in the equatorial ionosphere
> derived from CHAMP satellite magnetic field measurements, J. Atmos. Sol.-Terr.
> Phys., 72(4), 319-326,
> [doi:10.1016/j.jastp.2009.02.006](https://doi.org/10.1016/j.jastp.2009.02.006)
>
> EEJ version 2 additionally uses the EUVAC proxy and lunar local time. Coefficient
> files (CHAMP) are those distributed with Alken's `eej` / `eef` source packages.

#### `Empirical.run_ppeefm1` (Manoj & Maus PPEF — web service, or local RTEEF with `nativelypackaged_code=True`)

> Manoj, C., and S. Maus (2012), A real-time forecast service for the ionospheric
> equatorial zonal electric field, Space Weather, 10, S09002,
> [doi:10.1029/2012SW000825](https://doi.org/10.1029/2012SW000825)
>
> Manoj, C., S. Maus, H. Luehr, and P. Alken (2008), Penetration characteristics of
> the interplanetary electric field to the daytime equatorial ionosphere,
> J. Geophys. Res., 113, A12310,
> [doi:10.1029/2008JA013381](https://doi.org/10.1029/2008JA013381)

#### `Empirical.run_rocsat`

> Fejer, B. G., J. W. Jensen, and S.-Y. Su (2008),
> Quiet time equatorial F region vertical plasma drift model derived from
> ROCSAT-1 observations, J. Geophys. Res., 113, A05304,
> [doi:10.1029/2007JA012801](https://doi.org/10.1029/2007JA012801)

#### `Empirical.run_scherliessfejer`

> Fejer, B. G., and Scherliess, L. (1997), Empirical models of storm time
> equatorial zonal electric fields, J. Geophys. Res., 102(A11), 24047-24056,
> [doi:10.1029/97JA02164](https://doi.org/10.1029/97JA02164)
>
> Scherliess, L., and B. G. Fejer (1997), Storm time dependence of equatorial
> disturbance dynamo zonal electric fields, J. Geophys. Res., 102(A11), 24037–24046,
> [doi:10.1029/97JA02165](https://doi.org/10.1029/97JA02165)
>
> Scherliess, L., and Fejer, B. G. (1999), Radar and satellite global
> equatorial F region vertical drift model, J. Geophys. Res., 104(A4), 6829-6842,
> [doi:10.1029/1999JA900025](https://doi.org/10.1029/1999JA900025)

#### `Empirical.run_weimer05` — Weimer high-latitude electric potential and field-aligned current

> Weimer, D. R. (2005a), Improved ionospheric electrodynamic models and application to
> calculating Joule heating rates, J. Geophys. Res., 110, A05306,
> [doi:10.1029/2004JA010884](https://doi.org/10.1029/2004JA010884)
>
> Weimer, D. R. (2005b), Predicting surface geomagnetic variations using ionospheric
> electrodynamic models, J. Geophys. Res., 110, A12307,
> [doi:10.1029/2005JA011270](https://doi.org/10.1029/2005JA011270)

The bundled Fortran, translated by B. Foster (HAO/NCAR), is the spherical-cap-harmonic (SCHA) revision of the
2005a models that 2005b describes (degrees up to 12, orders up to 2, cap-size-dependent degrees, one shared
boundary). It gives the **electric potential** and the **field-aligned current** poleward of the model's
low-latitude boundary, driven by the IMF, the solar wind and the dipole tilt. `Empirical.run_weimer05` is a thin
wrapper over `models.weimer05` (same arguments, `obj.time` as the times).

> **Averaging default — please note.** Unless you pass the drivers yourself, `by`, `bz`, `vsw` and `nsw` come from
> the index store (NOAA OMNI) as the **mean of the previous 20 minutes of the 5-minute series**
> (`res='5min', avg=20`), which is how Weimer (2005b) drives the model — **not** as instantaneous values.
> Ask for those with `avg=None, res='1min'`; use `avg=45, lag=10` for the recipe the coefficients were fitted with
> (Weimer 2005a); drivers you pass in yourself are used as given, never averaged.

```python
import numpy as np
from datetime import datetime, timedelta
from mpyricalspace import Predictor, models

# one time step on a magnetic grid: a (mlat, mlt) cube, drivers given by hand
ds = models.weimer05(datetime(2024, 5, 11, 12), mlat=np.arange(50., 90.1, 1.), mlt=np.arange(0., 24., 1.),
                     by=0., bz=-5., vsw=450., nsw=9.)
ds.epot            # kV      (mlat, mlt)
ds.fac             # uA/m^2  (mlat, mlt), positive = downward, at 110 km

# through the facade: obj.time are the times; drivers from the index store, DEFAULT = the mean of the previous
# 20 min of the 5-min OMNI series (2005b)
obj = Predictor.Empirical(); obj.set_period(d0, dn, freq=timedelta(minutes=5))
ds = obj.run_weimer05(mlat, mlt)                                              # (time, mlat, mlt)

# the recipe the coefficients were derived with (2005a): the mean from 55 to 10 min before each step
ds = obj.run_weimer05(mlat, mlt, res='1min', avg=45, lag=10)
# instantaneous 1-min values instead of an average
ds = obj.run_weimer05(mlat, mlt, res='1min', avg=None)

# a satellite track: equal-length time / mlat / mlt arrays are aligned samples
ds = models.weimer05(times, track_mlat, track_mlt)
```

`examples/plot_weimer05_fig02.py` reproduces the paper's Figure 2 (nine polar maps of the potential) and checks
the extremes against the values printed in it; it runs as part of `tests/test_examples.py`.

Along a satellite track or over a geographic grid, `survey.run_track` / `survey.run_grid` include it
(`models=["weimer05"]`; the AACGM conversion uses `aacgmv2`, installed with the package — see
[Many models at once](#many-models-at-once--along-a-track-or-over-a-grid)).

- **Inputs.** `time`; `mlat`, `mlt` (below); `by`, `bz` [nT, GSM], `vsw` [km/s], `nsw` [cm⁻³] as a scalar or one
  value per time — omitted, they come from `DataManager` (`by bz vsw nsw`); `tilt` [deg] — omitted,
  `models.dipole_tilt(time)` computes it (IGRF-14 dipole, checked against astropy's Sun to 0.005°); `res`, `avg`,
  `lag` — how store-provided drivers are averaged: the mean of the `avg` minutes (of the `res` series) that end
  `lag` minutes before `time`. **Defaults `res='5min', avg=20, lag=0`**, i.e. the previous 20 min of 5-min data
  (see the note above); `avg=None` gives the value at `time`.
- **Coordinates.** The model was built in **AACGM** (altitude-adjusted corrected geomagnetic) latitude and
  magnetic local time (2005a, section 2 and appendix A); converting geographic → AACGM is up to you when you call it directly (`aacgmv2`,
  an AACGM implementation installed with the package, does it; the survey uses it automatically). The paper notes that modified-apex and AACGM
  coordinates differ very little at high latitudes, so apex-type latitudes such as the quasi-dipole `mlat` from
  `run_hltwim` should be close at these latitudes — that closeness has not been quantified here.
  `fac` is a density on the magnetic grid: to map it onto a geographic grid it must be compensated for the
  non-uniform mapping between the two systems (Richmond, 1995, as the papers note).
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

##### What the papers say, and where it lives in the code

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

##### Validation against the paper's figures

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

##### Build notes

The model is built by meson together with the rest of the package (`src/meson.build`). Its compiler flags are
requirements of the model, not tuning: `-fdefault-real-8` (default `real` and unsuffixed literals need 8 bytes) and
`-fno-automatic` (locals such as `skip` are read before they are set and rely on zero-initialised static storage).
Only the driver `weimer05_batch` is exposed to python; the model's own routines are compiled in and called from it. The `.dat` files are read from `_data/weimer05/` (an installed
package) or `src/weimer05/` (editable install).

##### Credits and terms

The Fortran (`src/weimer05/`) is a Fortran 90 translation, by B. Foster (HAO/NCAR, 2008), of D. Weimer's IDL model. It ships with its own `LICENSE`, which is **separate from this package's MIT license and
applies to those files** (the source, `weimer05_batch.f90` aside, and the four `.dat` coefficient files). It is
installed next to the data (`_data/weimer05/LICENSE`) and must stay with any redistribution. Its full text:

> This implementation of the Weimer model is provided for
> non-commercial scientific research and educational purposes
> only. Redistribution and modification for those purposes are
> permitted, provided that appropriate attribution to the original
> author and publications is retained. Commercial use,
> redistribution for commercial purposes, or incorporation into
> commercial products requires separate permission from the
> copyright holder. The software is provided ‘as is,’ without
> warranty of any kind.

In practice: use it for research and teaching, cite Weimer (2005a, 2005b) and credit the Fortran translation
(B. Foster, HAO/NCAR), and ask for permission before any commercial use or redistribution.
