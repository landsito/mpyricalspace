'''
Run several empirical models over a satellite track or a coordinate grid and
merge what applies into one xr.Dataset.

    from mpyricalspace import survey

    ds = survey.run_track(times, lats, lons, alts_km)          # aligned samples (an ephemeris)
    ds = survey.run_grid(times=..., lats=..., lons=..., alts=...)   # a meshgrid

    survey.run_track(..., models=["iri", "msis"])              # models=None (default) runs them all
    survey.model_domains()                                     # the allowed names + where each applies

Each model is evaluated only where it is defined; everything else is NaN.
``ds.attrs['models']`` lists what ran, ``ds.attrs['skipped']`` what did not and
why.  Variables carry the model name (``igrf_B``, ``hwm_u``, ``iri_NMF2``,
``eej``, ``eef``, ``sf_qvdrift``, ``manoj_ppef`` ...).

The global models (IGRF, MSIS, HWM, IRI) are evaluated at every sample.  The
equatorial-electrodynamics models describe a quantity at the magnetic equator as
a function of longitude and local time -- NOT of the sample's latitude -- so they
are computed only within a few degrees of the dip equator (IGRF dip latitude, at
most ``equator_deg``) and only inside their altitude regime:

    sf     Scherliess-Fejer *quiet* vertical drift   F region  (~200-900 km), |mlat| <= 2.5
    rocsat ROCSAT-1 quiet vertical drift             F region,                |mlat| <= 3
    eej    Alken equatorial electrojet current       E region  (~90-130 km),  |mlat| <= 5
    eef    Alken equatorial zonal electric field     E-F       (~90-1000 km), |mlat| <= 10
    manoj  Manoj-Maus RTEEF prompt-penetration E     F region  (~200-900 km), |mlat| <= 5;
           local RTEEF, so ACE-limited to one calendar year in 2001-2007

Not here: the Scherliess-Fejer *storm* decomposition (prompt / dynamo) is a
Jicamarca model -- use ``run_scherliessfejer`` at a fixed longitude for that.
``jvdm1`` is a Peruvian-sector model -- use its own ``run_*``.
``sf`` and ``manoj`` are two different prompt-penetration models (climatological
storm statistics giving a drift [m/s] vs a solar-wind transfer function giving
the equatorial E-field [mV/m]); pick one.
'''
import warnings

import numpy as np
import xarray as xr

from mpyricalspace import models as _m

_t2dt = _m._t2dt
_DEFAULT_ALT_KM = 300.0
_JRO_LON = _m.JRO[1]


# --- per-model track adapters ----------------------------------------------
# each takes flat (t, lat, lon, alt) arrays of equal length + an options dict
# and returns the model's xr.Dataset on the 'time' axis.

def _a_igrf(t, lat, lon, alt, o):
    return _m.igrf(t, list(lat), list(lon), list(alt), version=o.get('igrf_version', 14))


def _a_msis(t, lat, lon, alt, o):
    return _m.msis(t, list(lat), list(lon), list(alt), rho_m3=o.get('rho_m3', True),
                   nativelypackaged_indices=o.get('nativelypackaged_indices', False))


def _a_hwm(t, lat, lon, alt, o):
    return _m.hwm(t, list(lat), list(lon), list(alt), version=o.get('hwm_version', 2014))


def _a_iri(t, lat, lon, alt, o):
    return _m.iri(t, list(lat), list(lon), list(alt), version=o.get('iri_version', 2026),
                  compute_Ni=o.get('iri_ions', False),
                  nativelypackaged_indices=o.get('nativelypackaged_indices', False))


def _a_hltwim(t, lat, lon, alt, o):
    return _m.hltwim(t, list(lat), list(lon))


def _a_rocsat(t, lat, lon, alt, o):
    return _m.rocsat_drift(t, lon=list(lon))


def _a_eej(t, lat, lon, alt, o):
    return _m.eej(t, lon=list(lon), version=2, model=o.get('eej_model', 'champ'))


def _a_eef(t, lat, lon, alt, o):
    return _m.eef(t, lon=list(lon))


def _a_sf(t, lat, lon, alt, o):
    'Scherliess-Fejer quiet drift only -- the storm decomposition is Jicamarca-specific'
    return _m.scherliess_fejer(t, lon=list(lon))[['qvdrift']].rename(qvdrift='sf_qvdrift')


def _a_manoj(t, lat, lon, alt, o):
    '''Manoj-Maus RTEEF prompt-penetration equatorial E-field [mV/m], local model.
       ppef is the ACE IEF Ey run through the TF.COF filter -- a single
       solar-wind-driven series, ~uniform in longitude near the dip equator, so it
       is just resampled onto the track times.  For it as a vertical drift at a
       fixed longitude use run_ppeefm1.'''
    dts = _t2dt(t)
    y0, y1 = dts.min().year, dts.max().year
    if not (y0 == y1 and y0 in _m._RTEEF_YEARS):
        raise RuntimeError("local RTEEF has ACE data for 2001-2007 only, one calendar year")
    d0 = dts.min()
    raw = _m._rteef(d0, dts.max(), _JRO_LON,                      # [dt, qef, tef, ppef]
                    nativelypackaged_indices=o.get('nativelypackaged_indices', False))
    dhr = np.array([(x - d0).total_seconds() / 60. for x in raw[:, 0]])
    hr  = np.array([(x - d0).total_seconds() / 60. for x in dts])
    return xr.Dataset({'manoj_ppef': ('time', np.interp(hr, dhr, raw[:, 3].astype(float)))},
                      coords={'time': np.asarray(dts, dtype='datetime64[s]')})


# name -> dict(kind, adapter, alt_km window, magnetic-latitude cap [deg], needs an altitude?)
#   kind 'global'   : evaluated at every sample; masked by the altitude window
#        'high_lat' : the model itself returns NaN for |QD magnetic latitude| < 40
#        'equator'  : masked to |IGRF dip latitude| <= min(equator_deg, maglat) AND the altitude window
_REGISTRY = {
    'igrf':   dict(kind='global',   adapter=_a_igrf,   alt_km=None,          maglat=None, needs_alt=True),
    'msis':   dict(kind='global',   adapter=_a_msis,   alt_km=(0., 1000.),   maglat=None, needs_alt=True),
    'hwm':    dict(kind='global',   adapter=_a_hwm,    alt_km=(0., 500.),    maglat=None, needs_alt=True),
    'iri':    dict(kind='global',   adapter=_a_iri,    alt_km=(60., 2000.),  maglat=None, needs_alt=True),
    'hltwim': dict(kind='high_lat', adapter=_a_hltwim, alt_km=None,          maglat=None, needs_alt=False),
    'sf':     dict(kind='equator',  adapter=_a_sf,     alt_km=(200., 900.),  maglat=2.5,  needs_alt=False),
    'rocsat': dict(kind='equator',  adapter=_a_rocsat, alt_km=(200., 900.),  maglat=3.0,  needs_alt=False),
    'eej':    dict(kind='equator',  adapter=_a_eej,    alt_km=(90., 130.),   maglat=5.0,  needs_alt=False),
    'eef':    dict(kind='equator',  adapter=_a_eef,    alt_km=(90., 1000.),  maglat=10.0, needs_alt=False),
    'manoj':  dict(kind='equator',  adapter=_a_manoj,  alt_km=(200., 900.),  maglat=5.0,  needs_alt=False),
}

_FULL_NAME = {'sf': 'scherliess_fejer (quiet)', 'manoj': 'manoj_maus (local RTEEF)',
              'eej': 'eej', 'eef': 'eef', 'rocsat': 'rocsat_drift'}


def model_domains():
    '''``{model: {kind, altitude_km, maglat_deg, model, note}}`` -- where each model
       in run_track / run_grid applies.'''
    note = {
        'global':   'evaluated at every sample',
        'high_lat': 'auroral / polar cap only (|QD magnetic latitude| > 40 deg); NaN elsewhere',
        'equator':  'a magnetic-equator quantity vs longitude; masked to the dip equator and its altitude regime',
    }
    return {name: {'kind': r['kind'], 'altitude_km': r['alt_km'], 'maglat_deg': r['maglat'],
                   'model': _FULL_NAME.get(name, name), 'note': note[r['kind']]}
            for name, r in _REGISTRY.items()}


def _prefix(name, ds):
    'prefix a model result with <name>_ unless the variable already carries the name'
    out = {}
    for k, v in ds.data_vars.items():
        key = k if (k == name or k.startswith(name + '_') or k.endswith('_' + name)) else name + '_' + k
        out[key] = v
    return out


def run_track(times, lats, lons, alts=None, models=None, equator_deg=20.0, **opts):
    '''
    Evaluate every applicable empirical model along a trajectory; merge into one
    xr.Dataset on a shared ``time`` axis (dim ``time``, one row per sample).

    times, lats, lons : 1-D, equal length -- the samples (a satellite track).
    alts   : altitudes [km]; a scalar broadcasts.  Omitted -> 300 km (with a
             warning) for the models that need it.
    models : which models to run.  ``None`` (the default) runs every applicable
             model; otherwise pass a list of names, e.g. ``["iri", "msis"]`` or
             ``["sf", "manoj", "eef"]``.  ``list(survey.model_domains())`` is the
             list of allowed names; ``survey.model_domains()`` also gives each
             model's kind, altitude window and dip-latitude cap.  A requested
             name that is not a run_track model lands in ``ds.attrs['skipped']``.
    equator_deg : upper bound on how far from the dip equator the equatorial
             models are still evaluated; each model also has a tighter physical
             cap (see model_domains()).
    **opts : igrf_version / hwm_version / iri_version / eej_model / iri_ions / rho_m3 /
             nativelypackaged_indices (iri / msis / manoj -- False: DataManager (default);
             True: each model's own bundled index file).

    Out-of-domain samples come back as NaN.  ``ds.attrs['models']`` lists what
    ran, ``ds.attrs['skipped']`` what did not and why.
    '''
    t = np.asarray(_t2dt(times), dtype='datetime64[s]')
    n = t.size
    lat = np.broadcast_to(np.asarray(lats, float), n).astype(float)
    lon = np.broadcast_to(np.asarray(lons, float), n).astype(float)
    alt = np.full(n, np.nan) if alts is None else np.broadcast_to(np.asarray(alts, float), n).astype(float)

    want = list(_REGISTRY) if models is None else [str(m) for m in models]
    need_alt = any(_REGISTRY.get(m, {}).get('needs_alt') for m in want)
    if need_alt and not np.isfinite(alt).any():
        warnings.warn("run_track: no altitude given; using %g km for the global models" % _DEFAULT_ALT_KM)
        alt = np.full(n, _DEFAULT_ALT_KM)
    alt_or_default = np.where(np.isfinite(alt), alt, _DEFAULT_ALT_KM)

    need_dip = any(_REGISTRY.get(m, {}).get('kind') == 'equator' for m in want)
    dip = None
    if need_dip:
        try:
            dip = np.asarray(_m.igrf(t, list(lat), list(lon), list(alt_or_default)).dip.data, float)
        except Exception as e:                                     # pragma: no cover
            warnings.warn("run_track: could not compute dip latitude (%s); equatorial models skipped" % e)

    ds = xr.Dataset(coords={'time': t, 'lat': ('time', lat), 'lon': ('time', lon), 'alt': ('time', alt)})
    if dip is not None:
        ds = ds.assign_coords(dip_lat=('time', dip))

    ran, skipped = [], {}
    for name in want:
        if name not in _REGISTRY:
            skipped[name] = 'not a run_track model (see model_domains())'
            continue
        r = _REGISTRY[name]
        if r['kind'] == 'equator' and dip is None:
            skipped[name] = 'dip latitude unavailable'
            continue
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', message='.*outside the .* local-time window.*')
                res = r['adapter'](t, lat, lon, alt_or_default, opts)
        except Exception as e:
            skipped[name] = '%s: %s' % (type(e).__name__, e)
            continue

        mask = np.ones(n, bool)
        if r['alt_km'] is not None:
            mask &= (alt >= r['alt_km'][0]) & (alt <= r['alt_km'][1])
        if r['kind'] == 'equator':
            mask &= np.abs(dip) <= min(equator_deg, r['maglat'])

        for k, v in _prefix(name, res).items():
            arr = np.asarray(v.data, float)
            if arr.shape[:1] != (n,):
                ds[k] = (v.dims, arr)
                continue
            arr = arr.copy()
            arr[~mask] = np.nan
            ds[k] = (('time',) + tuple('%s_d%d' % (name, i) for i in range(1, arr.ndim)), arr)
        ran.append(name)

    ds.attrs.update(models=ran, skipped=skipped, equator_deg=float(equator_deg))
    return ds


def run_grid(times=None, lats=None, lons=None, alts=None, models=None, equator_deg=20.0, **opts):
    '''
    Evaluate every applicable model on the meshgrid of ``times x lats x lons x alts``;
    return one xr.Dataset with those dims (singleton dims squeezed).

    Same knobs as run_track -- in particular ``models=None`` runs every
    applicable model and a list (``models=["iri", "igrf"]``) restricts it; see
    ``survey.model_domains()`` for the names.  A dense global grid is a lot of
    points -- IGRF in particular loops in Python -- so keep the resolution modest.
    '''
    from datetime import datetime
    if times is None:
        times = [datetime(2020, 1, 1)]
    if lats is None:
        lats = np.arange(-88., 88.1, 4.)
    if lons is None:
        lons = np.arange(-180., 180., 15.)
    if alts is None:
        alts = [_DEFAULT_ALT_KM]

    tt = np.asarray(_t2dt(times), dtype='datetime64[s]')
    la = np.atleast_1d(np.asarray(lats, float))
    lo = np.atleast_1d(np.asarray(lons, float))
    al = np.atleast_1d(np.asarray(alts, float))
    T, LA, LO, AL = np.meshgrid(tt, la, lo, al, indexing='ij')
    shape = T.shape

    flat = run_track(T.ravel(), LA.ravel(), LO.ravel(), AL.ravel(),
                     models=models, equator_deg=equator_deg, **opts)

    out = xr.Dataset(coords={'time': tt, 'lat': la, 'lon': lo, 'alt': al})
    dims = ('time', 'lat', 'lon', 'alt')
    for k, v in flat.data_vars.items():
        if v.data.shape[:1] == (T.size,):
            out[k] = (dims, v.data.reshape(shape))
    if 'dip_lat' in flat.coords:
        out = out.assign_coords(dip_lat=(dims, flat['dip_lat'].data.reshape(shape)))
    out.attrs.update(flat.attrs)
    return out.squeeze()
