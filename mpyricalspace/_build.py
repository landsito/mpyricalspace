'''
Which compiled Fortran extensions actually built and load?

``mpyricalspace`` ships ~20 f2py extension modules (one or more per empirical
model). A partial or stale build leaves some of them missing, and because
``mpyricalspace.models`` imports them all at module top, the failure surfaces far
from its cause. ``build_info()`` probes each one in isolation.

    >>> import mpyricalspace
    >>> mpyricalspace.build_info()                 # {ext: None if OK else "<error>"}
    >>> print(mpyricalspace.build_report())        # a table
    $ python -m mpyricalspace doctor               # the same table
    $ python -m mpyricalspace doctor --json        # machine-readable (the VS Code helper reads this)
'''
import json
import importlib

# extension module -> (what it powers, its src/ subdirectory)
# kept in sync with src/meson.build
EXTENSIONS = {
    "rocsatf":          ("ROCSAT-1 vertical drift (Fejer et al. 2008)",        "rocsat"),
    "scherliessfejerf": ("Scherliess-Fejer vertical drift (quiet + storm)",    "scherliessfejer"),
    "jvdm1":            ("Alken JULIA vertical-drift model (2009)",            "jvdm1"),
    "eejm1":            ("Alken equatorial electrojet, v1",                    "eejm/eejm1"),
    "eejm2":            ("Alken equatorial electrojet, v2",                    "eejm/eejm2"),
    "eefm1":            ("Alken equatorial electric field",                   "eefm1"),
    "ppeefm1":          ("Manoj & Maus prompt-penetration EEF (RTEEF)",        "ppeefm1"),
    "hwm14f":           ("Horizontal Wind Model 2014",                        "hwm/hwm14"),
    "hwm07f":           ("Horizontal Wind Model 2007",                        "hwm/hwm07"),
    "hwm93f":           ("Horizontal Wind Model 1993",                        "hwm/hwm93"),
    "hltwimf":          ("HL-TWiM high-latitude wind (Dhadly et al. 2019)",    "hltwim"),
    "igrf14f":          ("IGRF-14 geomagnetic field",                         "igrf/igrf14"),
    "igrf13f":          ("IGRF-13 geomagnetic field",                         "igrf/igrf13"),
    "igrf12f":          ("IGRF-12 geomagnetic field",                         "igrf/igrf12"),
    "igrf11f":          ("IGRF-11 geomagnetic field",                         "igrf/igrf11"),
    "igrf10f":          ("IGRF-10 geomagnetic field",                         "igrf/igrf10"),
    "igrf09f":          ("IGRF-9 geomagnetic field",                          "igrf/igrf09"),
    "iri26f":           ("IRI-2026 ionosphere",                              "iri/iri26"),
    "iri20f":           ("IRI-2020 ionosphere",                              "iri/iri20"),
    "iri16f":           ("IRI-2016 ionosphere",                              "iri/iri16"),
    "iri12f":           ("IRI-2012 ionosphere",                              "iri/iri12"),
    "iri07f":           ("IRI-2007 ionosphere",                              "iri/iri07"),
    "iri01f":           ("IRI-2001 ionosphere",                              "iri/iri01"),
}


def build_info():
    '''Import every compiled extension on its own. Returns an ordered dict
    ``{ext: None}`` when it loads, ``{ext: "ErrorType: message"}`` when it does not.'''
    out = {}
    for ext in EXTENSIONS:
        try:
            importlib.import_module("mpyricalspace." + ext)
            out[ext] = None
        except Exception as e:                     # ImportError, or a dlopen OSError
            out[ext] = "%s: %s" % (type(e).__name__, e)
    return out


def build_report(info=None):
    '''A printable table of :func:`build_info`; ends with a one-line summary.'''
    info = build_info() if info is None else info
    w = max(len(e) for e in EXTENSIONS)
    lines = []
    for ext, (label, _sub) in EXTENSIONS.items():
        err = info.get(ext, "not probed")
        mark = "  ok  " if err is None else " FAIL "
        lines.append("  %-*s [%s] %s" % (w, ext, mark, label if err is None else err))
    ok = sum(v is None for v in info.values())
    lines.append("")
    lines.append("  %d/%d extensions OK" % (ok, len(EXTENSIONS)))
    return "\n".join(lines)


def build_json(info=None):
    '''JSON string: ``{"extensions": {ext: {ok, error, label, subdir}}, "ok": n, "total": m}``.'''
    info = build_info() if info is None else info
    exts = {ext: {"ok": info.get(ext) is None, "error": info.get(ext),
                  "label": label, "subdir": sub}
            for ext, (label, sub) in EXTENSIONS.items()}
    ok = sum(v["ok"] for v in exts.values())
    return json.dumps({"extensions": exts, "ok": ok, "total": len(exts)}, indent=2)
