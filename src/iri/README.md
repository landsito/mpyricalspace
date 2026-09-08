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

`sync_iri.py` runs automatically at build time (`meson setup`/`pip install`),
before compiling any `iriXX` extension, for all six versions.

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
