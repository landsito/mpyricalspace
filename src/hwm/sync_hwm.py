"""
sync_hwm.py -- downloads/updates hwm14 and hwm07's own files from
map.nrl.navy.mil. Runs automatically at build time (meson's
run_command(), see ../meson.build), before compiling either extension.
hwm93 has no NRL page (coefficients are DATA statements in the vendored
.f itself) and is not touched by this script.

For each file: compares Last-Modified/Content-Length (HEAD request)
against _manifest.json, downloads only if different. No network / site
down / blocked -> warning printed, existing vendored copy kept, build
never fails. NRL's server returns 403 for every *.dat request regardless
of whether the file exists (confirmed against readable directory
listings and a real HTML "not found" page reused for that response) --
dwm07b104i.dat/dwm07b_104i.dat and gd2qd.dat can therefore never be
synced this way; their committed copies (fetched once via
web.archive.org, byte-identical between hwm14 and hwm07's copies of the
same DWM07B/apex-grid files) are the only source and always kept as-is.

POST_PATCHES: compile/runtime fixes applied to a file right after
downloading it, so they survive the next resync. Empty for now -- both
hwm14.f90 and hwm07.01d.f90, as currently served by NRL, compile and run
unmodified against this package's own f2py wrapper (hwm14_batch.f90 /
hwm07_batch.f90).
"""
import json
import os
import sys
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST_PATH = os.path.join(HERE, "_manifest.json")
UA = "Mozilla/5.0 (compatible; mpyricalspace-hwm-sync/1.0)"
TIMEOUT = 15

# version dir -> (base URL, [filenames as published there]).
VERSIONS = {
    "hwm14": (
        "https://map.nrl.navy.mil/map/pub/nrl/HWM/HWM14/HWM14_ess224-sup-0002-supinfo",
        ["hwm14.f90", "checkhwm14.f90", "dwm07b104i.dat", "gd2qd.dat",
         "hwm123114.bin", "README.txt"],
    ),
    "hwm07": (
        "https://map.nrl.navy.mil/map/pub/nrl/HWM/HWM07",
        ["hwm07.01d.f90", "checkhwm07.f90", "dwm07b_104i.dat", "gd2qd.dat",
         "hwm071308e.dat", "readme.txt"],
    ),
}

# (version, filename) -> [(old substring, new substring), ...], applied in
# order, right after that file is freshly downloaded. See module docstring.
POST_PATCHES = {}


def _apply_post_patches(vdir, fn, text, out):
    patches = POST_PATCHES.get((vdir, fn))
    if not patches:
        return text
    for old, new in patches:
        if old in text:
            text = text.replace(old, new)
            print(f"sync_hwm: post-patch {vdir}/{fn}: {old!r} -> {new!r}", file=out)
    return text


def _http_head(url):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return r.headers.get("Last-Modified"), r.headers.get("Content-Length")


def _http_get(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return r.read()


def sync(versions=None, out=sys.stderr):
    """Sync `versions` (default: all of VERSIONS) against map.nrl.navy.mil.
    Returns True if anything changed on disk."""
    versions = versions or list(VERSIONS)
    manifest = {}
    if os.path.exists(MANIFEST_PATH):
        with open(MANIFEST_PATH) as f:
            manifest = json.load(f)

    changed = False
    for vdir in versions:
        base, files = VERSIONS[vdir]
        outdir = os.path.join(HERE, vdir)
        os.makedirs(outdir, exist_ok=True)
        for fn in files:
            key = f"{vdir}/{fn}"
            url = f"{base}/{fn}"
            local_path = os.path.join(outdir, fn)
            try:
                lm, cl = _http_head(url)
            except (urllib.error.URLError, OSError, TimeoutError) as e:
                print(f"sync_hwm: WARNING no se pudo verificar {key} ({e}) "
                      "-- se mantiene la copia local", file=out)
                continue

            if os.path.exists(local_path) and manifest.get(key) == [lm, cl]:
                continue  # unchanged since last sync

            try:
                data = _http_get(url)
            except (urllib.error.URLError, OSError, TimeoutError) as e:
                print(f"sync_hwm: WARNING no se pudo descargar {key} ({e}) "
                      "-- se mantiene la copia local", file=out)
                continue

            if fn.endswith((".f90", ".txt")):
                text = _apply_post_patches(vdir, fn, data.decode("latin-1"), out)
                data = text.encode("latin-1")

            with open(local_path, "wb") as f:
                f.write(data)
            manifest[key] = [lm, cl]
            changed = True
            print(f"sync_hwm: actualizado {key}", file=out)

    if changed:
        with open(MANIFEST_PATH, "w") as f:
            json.dump(manifest, f, indent=2, sort_keys=True)
            f.write("\n")
    return changed


if __name__ == "__main__":
    sync()
