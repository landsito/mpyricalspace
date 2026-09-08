# IGRF source tree

```
src/igrf/
  igrf09/ .. igrf14/   Each version's own igrf<N>.f (self-contained --
                       coefficients are DATA statements, no external files)
                       plus igrf<ver>f.pyf.
  sync_igrf.py         Downloads/updates igrf<N>.f from ngdc.noaa.gov.
  build_pyf.py         Regenerates a version's .pyf.
  _manifest.json       Tracks what sync_igrf.py last fetched.
```

Versions go back to `igrf09` (IGRF-9). Older releases (`igrf8` and earlier)
were never published as a standalone `.f` file, only as coefficient tables.

`sync_igrf.py` runs automatically at build time (`meson setup`/`pip install`),
before compiling any `igrfNN` extension. It only re-downloads a file if it
changed on ngdc.noaa.gov; no network just means the last synced copy is used.

## Two fixes applied on top of NOAA's own files

1. **`.pyf` intent annotations.** NOAA's `igrf<N>.f` doesn't declare
   `intent(in)`/`intent(out)` for `igrf<N>syn`'s 10 arguments, so a plain
   `f2py -h` scan produces the wrong calling convention. `build_pyf.py` adds
   the correct annotations to the *generated `.pyf`* -- the `.f` file itself
   is left untouched, matching NOAA exactly.

2. **`igrf9.f` / `igrf10.f` compile fix.** Both fail to build with a modern
   gfortran (`SIGN(1.1,X)` type mismatch in unrelated demo code, unrelated
   to `igrf<N>syn`). `sync_igrf.py` patches this automatically after
   downloading -- see `POST_PATCHES` in that file.
