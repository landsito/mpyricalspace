"""
sync_iri.py -- downloads/updates each version's .for source and dgrf/igrf/
mcsat .dat files from irimodel.org. Runs automatically at build time
(meson's run_command(), see ../meson.build), before compiling any iriXX
extension.

Does NOT touch iri/common/ (CCIR/URSI, apf107.dat/ig_rz.dat): CCIR/URSI are
static (unchanged since 2013); apf107.dat/ig_rz.dat are refreshed instead by
mpyricalspace.models.refresh_iri_indices() at runtime (ig_rz.dat always;
apf107.dat only when nativelypackaged_indices=True).

For each file: compares Last-Modified/Content-Length (HEAD request) against
_manifest.json, downloads only if different. All files for one version are
checked/re-fetched together in the same pass (irimodel.org republishes a
version's whole .for set together, same Last-Modified timestamp on every
file). No network / site down -> warning printed, existing vendored copy
kept, build never fails.

POST_PATCHES: compile/runtime fixes applied to a file right after
downloading it, so they survive the next resync (a one-time sed on the
checked-out file would get silently overwritten the next time
irimodel.org's own copy changes). Each version's irifun.for, as served
today, has fixed-size arrays too small for the index files this codebase
uses -- ionoindx/indrz (fed by ig_rz.dat, 853 monthly values needed) and
aap/af107 (fed by apf107.dat, 40000+ daily records needed). Sizes below are
each version's own starting point on irimodel.org; iri20/iri26's
ionoindx/indrz (1600) and iri07/iri01 (no aap/af107 at all) don't need a
patch for that particular array.
"""
import json
import os
import sys
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST_PATH = os.path.join(HERE, "_manifest.json")
BASE = "https://www.irimodel.org"
UA = "Mozilla/5.0 (compatible; mpyricalspace-iri-sync/1.0)"
TIMEOUT = 15

# irimodel.org serves this one file with a stray "(1)" suffix -- the plain
# name 404s.
REMOTE_NAME_OVERRIDES = {
    ("iri12", "igrf.for"): "igrf(1).for",
}

# version dir -> (site URL slug, [filenames as published on that page]).
# Deliberately excludes CCIR/URSI/apf107/ig_rz -- see module docstring.
VERSIONS = {
    "iri12": ("IRI-2012", [
        "cira.for", "igrf.for", "iridreg.for", "iriflip.for", "irifun.for",
        "irisub.for", "iritec.for", "iritest.for", "iriorbit.for", "iriorbitmax.for",
        "igrf2015.dat", "igrf2015s.dat",
    ] + [f"dgrf{y}.dat" for y in range(1945, 2011, 5)]),
    "iri16": ("IRI-2016", [
        "cira.for", "igrf.for", "iridreg.for", "iriflip.for", "irifun.for",
        "irisub.for", "iritec.for", "iritest.for", "irirtam.for", "irirtam-test.for",
        "igrf2020.dat", "igrf2020s.dat",
    ] + [f"dgrf{y}.dat" for y in range(1945, 2016, 5)]
      + [f"mcsat{m}.dat" for m in range(11, 23)]),
    "iri20": ("IRI-2020", [
        "cira.for", "igrf.for", "iridreg.for", "iriflip.for", "irifun.for",
        "irisub.for", "iritec.for", "iritest.for", "rocdrift.for",
        "irirtam.for", "irirtam-test.for",
        "igrf2025.dat", "igrf2025s.dat",
    ] + [f"dgrf{y}.dat" for y in range(1945, 2021, 5)]
      + [f"mcsat{m}.dat" for m in range(11, 23)]),
    "iri26": ("IRI-2026", [
        "cira.for", "igrf.for", "iridreg.for", "iriflip.for", "irifun.for",
        "irisub.for", "iritec.for", "iritest.for", "rocdrift.for",
        "igrf2025.dat", "igrf2025s.dat", "ibp_emp_coeffs.dat",
    ] + [f"dgrf{y}.dat" for y in range(1945, 2021, 5)]
      + [f"mcsat{m}.dat" for m in range(11, 23)]),
    "iri07": ("IRI-2007", [
        "cira.for", "igrf.for", "iridreg.for", "irifun.for", "irisub.for",
        "iritec.for", "iritest.for", "igrf2010.dat", "igrf2010s.dat",
    ] + [f"dgrf{y}.dat" for y in range(1945, 2006, 5)]),
    "iri01": ("IRI-2001", [
        "cira.for", "igrf.for", "iridreg.for", "irifun.for", "irisub.for",
        "iritec.for", "iritest.for", "igrf10.dat", "igrf10s.dat",
    ] + [f"dgrf{y}.dat" for y in
         ("45", "50", "55", "60", "65", "70", "75", "80", "85", "90", "95", "00", "05")]),
}

# (version, filename) -> [(old substring, new substring), ...], applied in order,
# right after that file is freshly downloaded. See module docstring.
POST_PATCHES = {
    # ionoindx/indrz: array used by tcon() to hold ig_rz.dat's data.
    # aig/arz: same COMMON /igrz/ slot, declared under different names by
    # read_ig_rz() -- must be bumped to the same size as ionoindx/indrz, or
    # iymst/iymend (also in that COMMON block) land at the wrong offset.
    # aap/af107: array used to hold apf107.dat's data.
    ("iri12", "irifun.for"): [
        ("ionoindx(806),indrz(806)", "ionoindx(2048),indrz(2048)"),
        ("aig(806),arz(806)", "aig(2048),arz(2048)"),
        ("aap(23000,9)", "aap(40000,9)"),
        ("af107(23000,3)", "af107(40000,3)"),
    ],
    ("iri16", "irifun.for"): [
        ("aig(806),arz(806)", "aig(2048),arz(2048)"),
        ("ionoindx(806),indrz(806)", "ionoindx(2048),indrz(2048)"),
        ("aap(27000,9)", "aap(40000,9)"),
        ("af107(27000,3)", "af107(40000,3)"),
    ],
    ("iri20", "irifun.for"): [
        # ionoindx/indrz and aig/arz are already 1600 here (853 needed).
        ("aap(27000,9)", "aap(40000,9)"),
        ("af107(27000,3)", "af107(40000,3)"),
    ],
    ("iri26", "irifun.for"): [
        ("aap(27000,9)", "aap(40000,9)"),
        ("af107(27000,3)", "af107(40000,3)"),
    ],
    ("iri07", "irifun.for"): [
        # No aap/af107 in this version -- uses a different F10.7 lookup.
        ("ionoindx(722),indrz(722)", "ionoindx(2048),indrz(2048)"),
    ],
    ("iri01", "irifun.for"): [
        ("ionoindx(722),indrz(722)", "ionoindx(2048),indrz(2048)"),
    ],
}


def _remote_url(vdir, fn):
    name = REMOTE_NAME_OVERRIDES.get((vdir, fn), fn)
    slug = VERSIONS[vdir][0]
    # percent-encode the handful of characters these overrides ever use
    name = name.replace("(", "%28").replace(")", "%29")
    return f"{BASE}/{slug}/{name}"


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
            print(f"sync_iri: post-patch {vdir}/{fn}: {old!r} -> {new!r}", file=out)
        # else: already patched, or the site fixed it upstream -- nothing to do.
    return text


def sync(versions=None, out=sys.stderr):
    """Sync `versions` (default: all of VERSIONS) against irimodel.org.
    Returns True if anything changed on disk."""
    versions = versions or list(VERSIONS)
    manifest = {}
    if os.path.exists(MANIFEST_PATH):
        with open(MANIFEST_PATH) as f:
            manifest = json.load(f)

    changed = False
    for vdir in versions:
        _slug, files = VERSIONS[vdir]
        outdir = os.path.join(HERE, vdir)
        os.makedirs(outdir, exist_ok=True)
        for fn in files:
            key = f"{vdir}/{fn}"
            url = _remote_url(vdir, fn)
            local_path = os.path.join(outdir, fn)
            try:
                lm, cl = _http_head(url)
            except (urllib.error.URLError, OSError, TimeoutError) as e:
                print(f"sync_iri: WARNING no se pudo verificar {key} ({e}) "
                      "-- se mantiene la copia local", file=out)
                continue

            if os.path.exists(local_path) and manifest.get(key) == [lm, cl]:
                continue  # unchanged since last sync

            try:
                data = _http_get(url)
            except (urllib.error.URLError, OSError, TimeoutError) as e:
                print(f"sync_iri: WARNING no se pudo descargar {key} ({e}) "
                      "-- se mantiene la copia local", file=out)
                continue

            if (vdir, fn) in POST_PATCHES:
                text = _apply_post_patches(vdir, fn, data.decode("latin-1"), out)
                data = text.encode("latin-1")

            with open(local_path, "wb") as f:
                f.write(data)
            manifest[key] = [lm, cl]
            changed = True
            print(f"sync_iri: actualizado {key}", file=out)

    if changed:
        with open(MANIFEST_PATH, "w") as f:
            json.dump(manifest, f, indent=2, sort_keys=True)
            f.write("\n")
    return changed


if __name__ == "__main__":
    sync()
