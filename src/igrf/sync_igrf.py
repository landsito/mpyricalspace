"""
sync_igrf.py -- downloads/updates each version's igrf<N>.f from
ngdc.noaa.gov. Runs automatically at build time (meson's run_command(), see
../meson.build), before compiling any igrfNN extension.

Each IGRF version is one self-contained file: coefficients are DATA
statements, no external data files, no files shared between versions.
Files live flat at https://www.ngdc.noaa.gov/IAGA/vmod/igrf<N>.f. Only
igrf14.f is linked from NOAA's current IGRF page; igrf9.f-igrf13.f are
still hosted there but unlinked. igrf8.f and earlier don't exist as a
standalone .f file (only raw coefficient tables) -- igrf09 is the oldest
version this repo supports.

For each file: compares Last-Modified/Content-Length (HEAD request) against
_manifest.json, downloads only if different. No network / site down ->
warning printed, existing vendored copy kept, build never fails.

POST_PATCHES: compile fixes applied to a file right after downloading it,
so they survive the next resync. igrf9.f/igrf10.f don't compile with a
modern gfortran (`SIGN(1.1,X)` type mismatch in unrelated demo code;
igrf11.f onward already has `SIGN(1.1D0,X)`).
"""
import json
import os
import sys
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST_PATH = os.path.join(HERE, "_manifest.json")
BASE = "https://www.ngdc.noaa.gov/IAGA/vmod"
UA = "Mozilla/5.0 (compatible; mpyricalspace-igrf-sync/1.0)"
TIMEOUT = 15

# version dir -> remote filename (as published at BASE/<filename>).
VERSIONS = {
    "igrf09": "igrf9.f",
    "igrf10": "igrf10.f",
    "igrf11": "igrf11.f",
    "igrf12": "igrf12.f",
    "igrf13": "igrf13.f",
    "igrf14": "igrf14.f",
}

# (version, filename) -> [(old substring, new substring), ...], applied in
# order, right after that file is freshly downloaded. See module docstring.
POST_PATCHES = {
    ("igrf09", "igrf9.f"): [("SIGN(1.1,X)", "SIGN(1.1D0,X)")],
    ("igrf10", "igrf10.f"): [("SIGN(1.1,X)", "SIGN(1.1D0,X)")],
}


def _http_head(url):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return r.headers.get("Last-Modified"), r.headers.get("Content-Length")


def _http_get(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return r.read()


def _apply_post_patches(vdir, fn, text, out):
    patches = POST_PATCHES.get((vdir, fn))
    if not patches:
        return text
    for old, new in patches:
        if old in text:
            text = text.replace(old, new)
            print(f"sync_igrf: post-patch {vdir}/{fn}: {old!r} -> {new!r}", file=out)
        # else: already patched, or the site fixed it upstream -- nothing to do.
    return text


def sync(versions=None, out=sys.stderr):
    """Sync `versions` (default: all of VERSIONS) against ngdc.noaa.gov.
    Returns True if anything changed on disk."""
    versions = versions or list(VERSIONS)
    manifest = {}
    if os.path.exists(MANIFEST_PATH):
        with open(MANIFEST_PATH) as f:
            manifest = json.load(f)

    changed = False
    for vdir in versions:
        fn = VERSIONS[vdir]
        outdir = os.path.join(HERE, vdir)
        os.makedirs(outdir, exist_ok=True)
        key = f"{vdir}/{fn}"
        url = f"{BASE}/{fn}"
        local_path = os.path.join(outdir, fn)
        try:
            lm, cl = _http_head(url)
        except (urllib.error.URLError, OSError, TimeoutError) as e:
            print(f"sync_igrf: WARNING no se pudo verificar {key} ({e}) "
                  "-- se mantiene la copia local", file=out)
            continue

        if os.path.exists(local_path) and manifest.get(key) == [lm, cl]:
            continue  # unchanged since last sync

        try:
            data = _http_get(url)
        except (urllib.error.URLError, OSError, TimeoutError) as e:
            print(f"sync_igrf: WARNING no se pudo descargar {key} ({e}) "
                  "-- se mantiene la copia local", file=out)
            continue

        text = _apply_post_patches(vdir, fn, data.decode("latin-1"), out)
        data = text.encode("latin-1")

        with open(local_path, "wb") as f:
            f.write(data)
        manifest[key] = [lm, cl]
        changed = True
        print(f"sync_igrf: actualizado {key}", file=out)

    if changed:
        with open(MANIFEST_PATH, "w") as f:
            json.dump(manifest, f, indent=2, sort_keys=True)
            f.write("\n")
    return changed


if __name__ == "__main__":
    sync()
