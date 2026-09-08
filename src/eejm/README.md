# EEJM source tree

```
src/eejm/
  eejm1/    eejm2/
            P. Alken's own source + champ/oer/sac_{mean,stddev}_coeffs
            (eejm1 has no oer/sac files -- CHAMP only), plus eejm<N>.c
            (our own CPython extension entry point) and Makefile (our
            own, builds that extension instead of the standalone
            eej_plot tool).
  sync_eejm.py      Downloads/updates eejm1/eejm2's own files from
                     geomag.colorado.edu.
  _manifest.json    Tracks what sync_eejm.py last fetched.
```

`sync_eejm.py` runs automatically at build time (`meson setup`/`pip
install`), before compiling `eejm1`/`eejm2`. Each version ships as one
archive (eejm1: a .zip, eejm2: a .tar.gz); a Last-Modified/Content-Length
HEAD check on the archive itself decides whether to re-download.

Not extracted: the archive's own Makefile, and eejm1's Windows
`eej_plot.exe`. eejm2's `gsl_multifit_ndlinear.h` is placed at
`ndlinear/gsl_multifit_ndlinear.h`, matching how this repo's other two
copies of P. Alken's ndlinear code (`src/eefm1/`, `src/ppeefm1/`) are
laid out.

Both eejm1/eej_basis.c and eejm2/eej_basis.c read GSL's bspline
workspace's `n` field directly; current GSL keeps that struct opaque,
so `sync_eejm.py`'s `POST_PATCHES` rewrites it to `gsl_bspline_ncoeffs()`
after every download, and likewise re-applies eejm2's ndlinear.c include
path fix.

`mpyricalspace.models._alken_ef()` chdir's into each version's own
folder (`_datadir()`/`fortran_cwd()`) before calling in, since the
compiled extension opens its coefficient files by plain name.
