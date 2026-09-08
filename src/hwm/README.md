# HWM source tree

```
src/hwm/
  hwm14/    hwm07/    hwm93/
            Each version's own .f90 source, .bin/.dat coefficient files,
            <ver>_batch.f90 (our own f2py driver) and <ver>f.pyf
            (f2py signature file). hwm14/ and hwm07/ also keep an
            original/ subfolder -- the locally-modified files this repo
            used before syncing against NRL, kept for reference, not
            built.
  sync_hwm.py       Downloads/updates hwm14/hwm07's own files from
                     map.nrl.navy.mil.
  _manifest.json    Tracks what sync_hwm.py last fetched.
```

`sync_hwm.py` runs automatically at build time (`meson setup`/`pip
install`), before compiling `hwm14f`/`hwm07f`. hwm93 has no NRL page
(coefficients are `DATA` statements in the vendored `.f` itself) and is
never touched by it.

NRL returns 403 for every `*.dat` request on this server, whether or not
the file exists (confirmed against the real directory listing and the
site's own "not found" page, reused verbatim for that response) --
`dwm07b104i.dat`/`dwm07b_104i.dat` and `gd2qd.dat` can never be synced
this way. Their committed copies (fetched once via web.archive.org) are
the only source; a sync attempt on them always warns and keeps the local
copy, same as any other unreachable file. `.f90`/`.bin`/`.txt` files sync
normally.

HWM07 and HWM14 both use NRL's DWM07B disturbance-wind model and its own
geodetic/quasi-dipole coordinate table -- `dwm07b104i.dat` (hwm14) and
`dwm07b_104i.dat` (hwm07) are byte-identical, likewise `gd2qd.dat`.

`mpyricalspace.models.hwm()`: HWM14 resolves its `.bin`/`.dat` files via
`$HWMPATH`; HWM07 opens its `.dat` files from the working directory, so
`_datadir()`/`fortran_cwd()` point it there; HWM93 has no data files.
