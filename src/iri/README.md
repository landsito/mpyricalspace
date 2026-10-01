# IRI source tree

```
src/iri/
  common/           CCIR/URSI (official COMMON_FILES page) + apf107.dat/ig_rz.dat
                     (seed copies from CHAIN/ECHAIM) -- shared by every version below.
  iri01/  iri07/  iri12/  iri16/  iri20/  iri26/
                     Each version's own .for source + dgrf/igrf/mcsat .dat files,
                     plus iri_batch.f90 (our own f2py driver, identical across
                     versions) and iri<ver>f.pyf (f2py signature file).
  sync_iri.py        Downloads/updates a version's .for/.dat set from irimodel.org.
  build_pyf.py       Regenerates a version's .pyf from its current .for set.
  _manifest.json     Tracks what sync_iri.py last fetched.
```

`sync_iri.py` runs only when `MPYRICALSPACE_SYNC_SOURCES=1` is set at build
time (`meson setup`/`pip install`) -- by default a build uses the committed
copy as is and never rewrites tracked files; run it by hand (`python
src/iri/sync_iri.py`) to check for updates, then review and commit the
result. It runs before compiling any `iriXX` extension, for all six
versions.

`common/` files are never duplicated into a version's own folder -- Fortran
opens files by plain name in the working directory, so
`mpyricalspace.models._iri_srcdirs()`/`_iri_workdir()` merge `common/` + the
version's own folder into one flat runtime directory at import time (a real
install already has them flattened together by `install_data()`).

## Fixes needed to build/run a version

Getting a fresh `.for`/`.dat` set to actually work (not just compile) takes
three fixes. All six versions need #2; only iri12/16 need the `aig`/`arz`
half of it.

1. **`.pyf` generation.** `irisub.for` has a malformed multi-block `COMMON`
   statement that breaks a plain `f2py -h` scan. Use `build_pyf.py`, not
   `f2py -h` directly -- see its own docstring.

2. **Fixed-size arrays too small for today's index files.** `irifun.for`
   declares `ionoindx`/`indrz` (fed by `ig_rz.dat`, needs 853 monthly
   values) and `aap`/`af107` (fed by `apf107.dat`, needs 40000+ daily
   records) with a hardcoded size, smaller than that in every version's
   current irimodel.org copy. `sync_iri.py`'s `POST_PATCHES` bumps these
   automatically after each download -- see that file for the exact sizes
   per version. `read_ig_rz()` and `tcon()` each declare their own copy of
   `ionoindx`/`indrz` under different names (`aig`/`arz` vs
   `ionoindx`/`indrz`, same `COMMON /igrz/` slot) -- both must be the same
   size, or `iymst`/`iymend` land at the wrong offset (same symptom as #3,
   different cause).

3. **`read_ig_rz()`/`readapf107()` are never called internally.**
   `irisub.for` calls `tcon()`, which reads `iymst`/`iymend` from `COMMON
   /igrz/` -- populated only by `read_ig_rz()`. Nothing in
   `irisub.for`/`irifun.for` calls it, and neither does IRI's own
   `iritest.for` (this is a gap in IRI-2012's own published bundle, not
   specific to this wrapper). Skipping the call raises no error: `iymst`/
   `iymend` stay 0, `tcon` rejects every date as "out of range", and every
   output is IRI's own `-1` "not computed" sentinel.
   `mpyricalspace.models.iri()` already calls both (gated on
   `hasattr(mod, 'read_ig_rz')`, since iri07/iri01 read both files inline
   instead and don't expose them as separate subroutines).

## IRI-2001 original data files

irimodel.org's `IRI-2001/` folder no longer matches its own `igrf.for` and
`00_iri2001_readme.txt`: in April 2013 its IGRF-2000 files were replaced by
IGRF-11 ones (`dgrf00.dat`, `dgrf05.dat`, `igrf10.dat`, `igrf10s.dat`), but
`igrf.for` still opens `igrf00.dat`/`igrf00s.dat`, and `ap.dat` (read by
`irifun.for`'s `APF()` for the storm model, `quiet=False`) is gone too. With
the site's files alone, IRI-2001 hits a Fortran `STOP` for any date from 1995
on and a runtime error in storm mode -- both end the host Python process.

The three originals are vendored in `iri01/` and are **not** managed by
`sync_iri.py` (not in its `VERSIONS` list), so a resync never touches them:

| file | content | source |
|---|---|---|
| `igrf00.dat` | IGRF-2000 main field, degree 10 | AGI mirror, 2005-10-27 |
| `igrf00s.dat` | IGRF-2000 secular variation (2000-2005), degree 8 | AGI mirror, 2005-10-27 |
| `ap.dat` | daily 3-h ap + F10.7, 1960-01-01 .. 2006-10-31 | AGI mirror; CRLF -> LF |

Source: <https://ftp.agi.com/pub/STKData/CentralBodies/Earth/Rf/IRI2001ModelData/>.
The other files there are byte-identical to ours (CCIR/URSI, `dgrf45`-`dgrf90`)
or superseded (`ig_rz.dat`: ours is newer and refreshed at runtime;
`dgrf95.dat`: AGI's is the provisional IGRF-1995, ours the definitive DGRF-1995
irimodel.org now serves). `ap.dat` was converted from CRLF to LF line endings
only: `APF()` reads it as direct-access records of `RECL=39` (38 characters +
LF), one per day, so CRLF would misalign every record after the first.
After 2006-10-31 `APF()` finds no record and the storm correction is skipped.

The IGRF-11 files irimodel.org now serves for IRI-2001 are not used by its
`igrf.for` and are not synced or installed.
