"""
sync_eejm.py -- downloads/updates eejm1/eejm2's own files from
geomag.colorado.edu. Runs automatically at build time (meson's
run_command(), see ../meson.build), before compiling either extension.

Each version ships as a single archive (eejm1: EEJM-1.0.zip, eejm2:
EEJM-2.0.tar.gz) with one top-level directory inside. Compares the
archive's own Last-Modified/Content-Length (HEAD request) against
_manifest.json, re-downloads/re-extracts only if different. No network
/ site down -> warning printed, existing vendored copy kept, build
never fails.

Not extracted: the archive's own Makefile (this repo has its own,
building the CPython extension meson installs instead of the
standalone eej_plot tool) and eejm1's eej_plot.exe (Windows binary, no
build use here).

FILE_REMAP: eejm2's gsl_multifit_ndlinear.h is placed at
ndlinear/gsl_multifit_ndlinear.h, not flat -- matching how this repo's
other two local copies of P. Alken's ndlinear code (src/eefm1/,
src/ppeefm1/) are laid out.

POST_PATCHES: compile fixes applied to a file right after extracting
it, so they survive the next resync. eejm1/eej_basis.c and
eejm2/eej_basis.c both access GSL's bspline workspace's `n` field
directly; current GSL (the struct is opaque now) no longer exposes it
-- use gsl_bspline_ncoeffs() instead. eejm2/ndlinear.c's #include is
patched to match the ndlinear/ remap above.
"""
import io
import json
import os
import sys
import tarfile
import urllib.error
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST_PATH = os.path.join(HERE, "_manifest.json")
UA = "Mozilla/5.0 (compatible; mpyricalspace-eejm-sync/1.0)"
TIMEOUT = 15

# version dir -> (archive URL, archive kind, top-level dir name inside it).
VERSIONS = {
    "eejm1": ("https://geomag.colorado.edu/images/EEJM2/EEJM-1.0.zip", "zip", "EEJM-1.0"),
    "eejm2": ("https://geomag.colorado.edu/images/EEJM2/EEJM-2.0.tar.gz", "tar", "EEJM-2.0"),
}

EXCLUDE = {"Makefile", "eej_plot.exe"}

# (version, archive filename) -> local relative path.
FILE_REMAP = {
    ("eejm2", "gsl_multifit_ndlinear.h"): "ndlinear/gsl_multifit_ndlinear.h",
}

# (version, local relative path) -> [(old substring, new substring), ...],
# applied in order, right after that file is freshly extracted. See module
# docstring.
POST_PATCHES = {
    ("eejm1", "eej_basis.c"): [
        ("w->B = gsl_vector_alloc(w->bspline_workspace_p->n);",
         "w->B = gsl_vector_alloc(gsl_bspline_ncoeffs(w->bspline_workspace_p));"),
        ("gsl_vector_view v = gsl_vector_view_array(y, w->bspline_workspace_p->n);",
         "gsl_vector_view v = gsl_vector_view_array(y, gsl_bspline_ncoeffs(w->bspline_workspace_p));"),
    ],
    ("eejm2", "eej_basis.c"): [
        ("w->B = gsl_vector_alloc(w->bspline_workspace_p->n);",
         "w->B = gsl_vector_alloc(gsl_bspline_ncoeffs(w->bspline_workspace_p));"),
    ],
    ("eejm2", "ndlinear.c"): [
        ('#include "gsl_multifit_ndlinear.h"', '#include "ndlinear/gsl_multifit_ndlinear.h"'),
    ],
}


def _http_head(url):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return r.headers.get("Last-Modified"), r.headers.get("Content-Length")


def _http_get(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return r.read()


def _members(data, kind, topdir):
    """Yield (archive-relative filename, bytes) for every regular file in
    the archive, with the top-level directory stripped."""
    prefix = topdir + "/"
    if kind == "zip":
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            for info in zf.infolist():
                if info.is_dir() or not info.filename.startswith(prefix):
                    continue
                yield info.filename[len(prefix):], zf.read(info)
    else:
        with tarfile.open(fileobj=io.BytesIO(data)) as tf:
            for member in tf.getmembers():
                if not member.isfile() or not member.name.startswith(prefix):
                    continue
                yield member.name[len(prefix):], tf.extractfile(member).read()


def _apply_post_patches(vdir, relpath, text, out):
    patches = POST_PATCHES.get((vdir, relpath))
    if not patches:
        return text
    for old, new in patches:
        if old in text:
            text = text.replace(old, new)
            print(f"sync_eejm: post-patch {vdir}/{relpath}: {old!r} -> {new!r}", file=out)
    return text


def sync(versions=None, out=sys.stderr):
    """Sync `versions` (default: all of VERSIONS) against geomag.colorado.edu.
    Returns True if anything changed on disk."""
    versions = versions or list(VERSIONS)
    manifest = {}
    if os.path.exists(MANIFEST_PATH):
        with open(MANIFEST_PATH) as f:
            manifest = json.load(f)

    changed = False
    for vdir in versions:
        url, kind, topdir = VERSIONS[vdir]
        outdir = os.path.join(HERE, vdir)
        try:
            lm, cl = _http_head(url)
        except (urllib.error.URLError, OSError, TimeoutError) as e:
            print(f"sync_eejm: WARNING no se pudo verificar {vdir} ({e}) "
                  "-- se mantiene la copia local", file=out)
            continue

        if os.path.isdir(outdir) and manifest.get(vdir) == [lm, cl]:
            continue  # unchanged since last sync

        try:
            data = _http_get(url)
        except (urllib.error.URLError, OSError, TimeoutError) as e:
            print(f"sync_eejm: WARNING no se pudo descargar {vdir} ({e}) "
                  "-- se mantiene la copia local", file=out)
            continue

        os.makedirs(outdir, exist_ok=True)
        for name, raw in _members(data, kind, topdir):
            if name in EXCLUDE:
                continue
            relpath = FILE_REMAP.get((vdir, name), name)
            text = _apply_post_patches(vdir, relpath, raw.decode("latin-1"), out)
            local_path = os.path.join(outdir, relpath)
            os.makedirs(os.path.dirname(local_path) or outdir, exist_ok=True)
            with open(local_path, "wb") as f:
                f.write(text.encode("latin-1"))

        manifest[vdir] = [lm, cl]
        changed = True
        print(f"sync_eejm: actualizado {vdir}", file=out)

    if changed:
        with open(MANIFEST_PATH, "w") as f:
            json.dump(manifest, f, indent=2, sort_keys=True)
            f.write("\n")
    return changed


if __name__ == "__main__":
    sync()
