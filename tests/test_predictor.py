'''
Network-free tests: every model is exercised with indices passed explicitly so
nothing hits the DataManager / the web.
'''
from datetime import datetime, timedelta
import os
import numpy as np
import pytest
import xarray as xr
from mpyricalspace.Predictor import Empirical

JRO = (-11.958256, -76.859012, 250)


@pytest.fixture
def obj():
    o = Empirical()
    o.set_period(datetime(2024, 5, 11, 18), datetime(2024, 5, 11, 22), timedelta(hours=1))
    o.set_location(*JRO)
    return o


def test_igrf_no_indices(obj):
    ds = obj.run_igrf(version=13)
    assert isinstance(ds, xr.Dataset)
    assert set(ds.data_vars) >= {"Bx", "By", "Bz", "B", "dip", "dec", "inc"}
    assert np.all(np.isfinite(ds.B.data)) and np.all(ds.B.data > 0)      # |B| ~ 2.2e-5 T at JRO
    assert 1e-5 < np.mean(ds.B.data) < 6e-5


@pytest.mark.parametrize("version", [14, 13])
def test_igrf_versions_agree(obj, version):
    '''IGRF-14/13/12 share the same synthesis routine (just one more coefficient
       epoch in -14); they should agree closely at a shared, well-covered date.'''
    ds = obj.run_igrf(version=version)
    assert np.all(np.isfinite(ds.B.data))
    assert 1e-5 < np.mean(ds.B.data) < 6e-5
    if version == 14:
        assert np.mean(ds.B.data) == pytest.approx(float(obj.run_igrf(version=13).B.data.mean()), rel=1e-3)


def test_igrf14_is_default(obj):
    assert np.allclose(obj.run_igrf().B.data, obj.run_igrf(version=14).B.data)


def test_igrf14_secular_variation_extrapolation():
    '''2027 falls in IGRF-14's 2025-2030 provisional SV block (IGRF-13 only covers to 2025).'''
    o = Empirical()
    o.set_time([datetime(2027, 1, 1)])
    ds = o.run_igrf(version=14, lats=-12., lons=-77., alts=250.)
    assert np.isfinite(float(ds.B.data)) and 1e-5 < float(ds.B.data) < 6e-5


def test_iri_no_indices(obj):
    '''nativelypackaged_indices=True: F10.7 from IRI's own bundled apf107.dat, no DataManager hit.'''
    ds = obj.run_iri(alts=300., lats=-12., lons=-77., version=2020, nativelypackaged_indices=True)
    assert np.all(ds.ne.data > 1e9)                                       # F-region Ne, m^-3
    assert np.all((ds.Te.data > 300) & (ds.Te.data < 6000))


@pytest.mark.parametrize("version", [2012, 2016, 2020, 2026])
def test_iri_datamanager_is_default(version):
    '''nativelypackaged_indices=False (default): F10.7/F10.7a come from DataManager, same as
       every other model -- not from IRI's own apf107.dat.'''
    o = Empirical()
    o.set_time([datetime(2013, 6, 15, 12)])
    native = o.run_iri(alts=300., lats=-12., lons=-77., version=version, nativelypackaged_indices=True)
    default = o.run_iri(alts=300., lats=-12., lons=-77., version=version)
    assert 1e10 < float(default.ne.data) < 1e13
    # both are real F10.7 sources for the same date -- not asserting a difference, just that
    # the default path runs standalone (DataManager, not the bundled file) and returns sane output
    assert 300 < float(native.Te.data) < 6000


@pytest.mark.parametrize("version", [2012, 2016, 2020, 2026])
def test_iri_versions(version):
    o = Empirical()
    o.set_time([datetime(2013, 6, 15, 12)])                               # in range for all four
    ds = o.run_iri(alts=300., lats=-12., lons=-77., version=version)
    assert 1e10 < float(ds.ne.data) < 1e13
    assert 300 < float(ds.Te.data) < 6000


@pytest.mark.parametrize("version", [2016, 2020, 2026])
def test_iri_recent_date(version):
    '''ig_rz.dat (Rz12/IG12, always file-read) auto-refreshes past the packaged table (was:
       -1 overflow of the hardcoded aig(806)/ionoindx(806) arrays) with the default,
       DataManager-backed F10.7; nativelypackaged_indices=True additionally exercises the
       apf107.dat auto-refresh for F10.7 itself.'''
    o = Empirical()
    o.set_time([datetime(2026, 6, 15, 12)])
    ds = o.run_iri(alts=300., lats=-12., lons=-77., version=version)
    assert 1e10 < float(ds.ne.data) < 1e13
    ds_native = o.run_iri(alts=300., lats=-12., lons=-77., version=version, nativelypackaged_indices=True)
    assert 1e10 < float(ds_native.ne.data) < 1e13


def test_iri_external_f107():
    '''explicit F107 / F107a override DataManager/apf107.dat either way (topside Ne + Te
       scale with it); Rz12 / IG12 stay file-read so NmF2 is unchanged.'''
    o = Empirical()
    o.set_time([datetime(2013, 6, 15, 12)])
    lo = o.run_iri(alts=600., lats=-12., lons=-77., F107=70., F107a=70.)
    hi = o.run_iri(alts=600., lats=-12., lons=-77., F107=210., F107a=205.)
    assert float(hi.ne.data) > float(lo.ne.data) * 1.15
    assert float(hi.Te.data) > float(lo.Te.data) + 100
    assert abs(float(hi.NMF2.data) - float(lo.NMF2.data)) / float(lo.NMF2.data) < 0.02


def test_iri_bad_version(obj):
    with pytest.raises(ValueError):
        obj.run_iri(alts=300., lats=0., lons=0., version=2007)


def test_iri_leaves_cwd_untouched(obj):
    here = os.getcwd()
    obj.run_iri(alts=300., lats=0., lons=0.)
    assert os.getcwd() == here


def test_model_datadirs_exist():
    from mpyricalspace import models
    for name in ('hwm14', 'hwm07', 'hltwim', 'eejm1', 'eejm2', 'eefm1', 'ppeefm1',
                 'iri12', 'iri16', 'iri20', 'iri26'):
        assert os.path.isdir(models._datadir(name))
    assert os.environ.get('HWMPATH') == models._datadir('hwm14')


def test_eej_reference_row():
    '''eejm1 matches P. Alken's standalone eej_plot: lon 90, LT 10.5, doy 79, F10.7 180.'''
    o = Empirical()
    o.set_time([datetime(2024, 3, 20, 12)])
    ds = o.run_eej(version=1, lon=90., flux=180., slts=[10.5], doys=[79])
    assert float(ds.eej) == pytest.approx(0.155647, abs=1e-4)
    assert float(ds.eej_sigma) == pytest.approx(0.029475, abs=1e-4)


@pytest.mark.parametrize("version", [1, 2])
def test_eej_diurnal(version):
    o = Empirical()
    o.set_time([datetime(2024, 3, 20, 12)])
    j = np.array([float(o.run_eej(version=version, lon=-77., flux=150., doys=[80], slts=[t]).eej)
                  for t in (7., 11., 17.)])
    assert j[1] > j[0] and j[1] > j[2]                    # electrojet peaks near local noon
    assert np.all(np.abs(j) < 1.0)


def test_eej_outside_lt_window_is_nan():
    o = Empirical()
    o.set_time([datetime(2024, 3, 20, 12)])
    ds = o.run_eej(version=2, lon=0., flux=120., doys=[80], slts=[2.0], lunars=[12.0])
    assert np.isnan(float(ds.eej))


@pytest.mark.parametrize("model", ["champ", "oersted", "sac-c"])
def test_eej_v2_satellite_models(model):
    o = Empirical()
    o.set_time([datetime(2005, 1, 1)])
    lons = np.arange(-180.0, 181.0, 30.0)
    j = np.atleast_1d(o.run_eej(version=2, model=model, lon=lons, flux=150.,
                                slts=[10.5], doys=[80], lunars=12.4).eej.data)
    assert j.shape == lons.shape and np.all(np.abs(j) < 0.5)
    if model == "sac-c":                                  # SAC-C: J is local-time independent
        j2 = float(o.run_eej(version=2, model="sac-c", lon=0., flux=150.,
                             slts=[15.0], doys=[80], lunars=12.4).eej)
        j1 = float(o.run_eej(version=2, model="sac-c", lon=0., flux=150.,
                             slts=[10.5], doys=[80], lunars=12.4).eej)
        assert abs(j1 - j2) < 1e-9


def test_eej_rejects_bad_model():
    o = Empirical()
    o.set_time([datetime(2005, 1, 1)])
    with pytest.raises(ValueError):
        o.run_eej(version=2, model="viking", lon=0., flux=150., slts=[10.5], doys=[80])
    with pytest.raises(ValueError):
        o.run_eej(version=1, model="sac-c", lon=0., flux=150., slts=[10.5], doys=[80])


def test_eef_runs():
    o = Empirical()
    o.set_time([datetime(2024, 3, 20, 12)])
    ds = o.run_eef(lon=90., flux=100., slts=[10.5], doys=[266], lunars=[12.5])
    assert set(ds.data_vars) == {"eef", "eef_sigma"}
    assert 0.05 < float(ds.eef) < 5.0 and float(ds.eef_sigma) > 0


def test_fortran_batch_drivers_present():
    from mpyricalspace import (rocsatf, hwm14f, hwm07f, hwm93f, hltwimf,
                               iri26f, iri20f, iri16f, iri12f)
    assert hasattr(rocsatf, 'getverticaldrift_batch')
    assert hasattr(hwm14f, 'hwm14_batch')
    assert hasattr(hwm07f, 'hwm07_batch') and hasattr(hwm93f, 'gws5_batch')
    assert hasattr(hltwimf, 'hltwim_batch')
    assert all(hasattr(m, 'iri_batch') for m in (iri26f, iri20f, iri16f, iri12f))


def test_build_info_all_extensions_load():
    import mpyricalspace
    info = mpyricalspace.build_info()
    assert len(info) == len(mpyricalspace._build.EXTENSIONS)
    broken = {k: v for k, v in info.items() if v is not None}
    assert not broken, "extensions failed to import: %s" % broken
    n = len(info)
    assert "%d/%d extensions OK" % (n, n) in mpyricalspace.build_report(info)


def test_vscode_ext_vsix_builds(tmp_path):
    import zipfile
    from mpyricalspace import _vscode
    vsix = _vscode.build_vsix(str(tmp_path))
    with zipfile.ZipFile(vsix) as z:
        names = z.namelist()
    assert "extension.vsixmanifest" in names and "[Content_Types].xml" in names
    assert "extension/package.json" in names and "extension/extension.js" in names


def test_hltwim_reference_row():
    '''close to the checkhltwim.f90 reference (day 10, ut 21, kp 3, 65N 215E).
       ~0.3% off it: our initalf() bounds fix changes the coordinate transform
       slightly vs NRL's pasted numbers, which came from the unpatched (out-of-
       bounds) build.'''
    o = Empirical()
    o.set_time([datetime(2020, 1, 1)])
    ds = o.run_hltwim(kp=3., lat=65., lon=215., ut=21., doy=10)
    assert float(ds.v) == pytest.approx(45.492, abs=0.5)       # geographic meridional
    assert float(ds.u) == pytest.approx(-53.303, abs=0.5)      # geographic zonal
    assert float(ds.mv) == pytest.approx(20.675, abs=0.5)      # QD meridional
    assert float(ds.mu) == pytest.approx(-59.953, abs=0.5)     # QD zonal
    assert 55 < float(ds.mlat) < 75 and 0 <= float(ds.mlt) < 24  # 65N 215E ~ auroral, evening


def test_hltwim_low_latitude_is_nan(obj):
    ds = obj.run_hltwim(kp=3., lat=5., lon=0.)
    assert np.all(np.isnan(ds.u.data)) and np.all(np.isnan(ds.v.data))
    assert np.all(np.isfinite(ds.mlat.data)) and np.all(np.abs(ds.mlat.data) < 40)  # QD lat still reported


def test_hltwim_time_series(obj):
    ds = obj.run_hltwim(kp=4., lat=72., lon=-100.)
    assert ds.u.data.shape == obj.time.shape
    assert np.all(np.abs(ds.u.data) < 2000) and np.all(np.abs(ds.v.data) < 2000)


def test_hltwim_kp_dependence():
    '''storm winds exceed quiet winds over the polar cap -- the Kp=5 vs Kp=2
       contrast of Fig. 12 in Dhadly et al. (2019).'''
    from mpyricalspace import models
    glat, glon = np.meshgrid(np.arange(60., 89., 3.), np.arange(0., 360., 20.), indexing='ij')
    t = [datetime(2016, 3, 20, 8)]
    lo = models.hltwim(t, glat.ravel(), glon.ravel(), kp=2., ut=8., doy=80)
    hi = models.hltwim(t, glat.ravel(), glon.ravel(), kp=5., ut=8., doy=80)
    spd = lambda d: np.hypot(d.mu.data, d.mv.data)
    assert np.nanmean(spd(hi)) > 1.3 * np.nanmean(spd(lo))
    assert np.nanmax(np.abs(hi.mlat.data)) > 80                     # grid reaches the pole


def test_rocsat_explicit(obj):
    ds = obj.run_rocsat(f107s=np.full(obj.time.size, 150.))
    assert ds.qvdrift_rocsat.data.shape == obj.time.shape
    assert np.all(np.abs(ds.qvdrift_rocsat.data) < 200)


def test_rocsat_lon_doy_slt_grid(obj):
    lons, doys, slts = np.arange(-150., 151., 30.), np.arange(80, 111, 10), np.arange(0., 24., 1.)
    ds = obj.run_rocsat(f107s=150., lons=lons, doys=doys, slts=slts)
    assert ds.qvdrift_rocsat.dims == ("lon", "doy", "slt")
    assert ds.qvdrift_rocsat.shape == (lons.size, doys.size, slts.size)
    assert np.all(np.abs(ds.qvdrift_rocsat.data) < 200)


def test_hwm_explicit(obj):
    ds = obj.run_hwm(ap=5., lon=JRO[1], lat=JRO[0], alt=250., version=2014)
    assert np.all(np.abs(ds.u.data) < 800) and np.all(np.abs(ds.v.data) < 800)


@pytest.mark.parametrize("version", [2014, 2007, 1993])
def test_hwm_versions(obj, version):
    ds = obj.run_hwm(ap=5., lon=JRO[1], lat=JRO[0], alt=250., version=version)
    assert ds.u.data.shape == obj.time.shape
    assert np.all(np.abs(ds.u.data) < 1000) and np.all(np.abs(ds.v.data) < 1000)
    assert len(np.unique(np.round(ds.u.data, 3))) > 1        # varies with time, not a broadcast blob


def test_hwm_bad_version(obj):
    with pytest.raises(ValueError):
        obj.run_hwm(version=2020)


def test_msis_explicit(obj):
    ds = obj.run_msis(f107s=150., ap=5., lon=JRO[1], lat=JRO[0], alt=400.)
    assert np.all(ds.rho.data > 0) and np.all(ds.T.data > 400)


def test_jvdm1_explicit():
    o = Empirical()
    o.set_period(datetime(2024, 3, 20), datetime(2024, 3, 20, 4), timedelta(hours=1))
    ds = o.run_jvdm1(f107=100., f107a=100., slt=np.arange(9, 16), doy=80)
    assert "qvdrift_150km" in ds


def test_jvdm1_scalar_over_time(obj):
    # scalar f107/f107a gridded against the time axis (used to raise in the old Predictor)
    ds = obj.run_jvdm1(f107=100., f107a=100.)
    assert ds.qvdrift_150km.data.shape == obj.time.shape


def test_manoj_maus_with_explicit_data(monkeypatch):
    from mpyricalspace import models
    d0 = datetime(2024, 5, 11, 0)
    tt = np.array([d0 + timedelta(hours=i) for i in range(4)], dtype='datetime64[s]')
    monkeypatch.setattr(models, 'scherliess_fejer', lambda *a, **k: xr.Dataset(
        {'qvdrift': ('time', np.full(4, 20.))}, coords={'time': tt.astype('datetime64[s]').astype(datetime)}))
    web = np.array([[d0 + timedelta(hours=i), 0.5, 0.6, 0.1] for i in range(5)], dtype=object)
    ds = models.manoj_maus(tt, d0, d0 + timedelta(hours=3), timedelta(hours=1), data=web)
    assert set(ds.data_vars) == {"mV_to_ms", "vdrift_manoj", "prompt_manoj", "qvdrift"}
    assert ds.vdrift_manoj.data.shape == (4,)


def test_manoj_maus_synthetic_ief_transfer_function():
    """run_ppeefm1(ief=...) pushes a synthetic IEF Ey step through the TF.COF filter."""
    obj = Empirical()
    d0 = datetime(2004, 6, 1, 12)
    obj.set_period(d0, d0 + timedelta(hours=6), freq=timedelta(minutes=5))
    n = obj.time.size
    ief = np.zeros(n + 12)
    ief[12 + 12:] = 1.0                                  # 1 mV/m step, 1 h after the start
    ds = obj.run_ppeefm1(nativelypackaged_code=True, ief=ief)
    assert set(ds.data_vars) == {"qef", "ppef", "tef"}
    assert ds.ppef.size == n
    assert np.allclose(ds.qef.data, 0.0)                 # climatology skipped in synthetic mode
    pp = ds.ppef.data
    assert pp[0] == 0.0 and 0.03 < pp.max() < 0.06       # transient overshoot ~0.043 mV/m
    assert abs(pp[-1]) < 0.01                            # decays toward the small DC gain (shielding)
    with pytest.raises(ValueError):                      # wrong length is rejected
        obj.run_ppeefm1(nativelypackaged_code=True, ief=np.zeros(n))


def _sf_indices(d0, dn, ae=250., f107=140., step_min=15):
    '''synthetic [datetime, AE, F10.7] covering d0-30h .. dn for scherliess_fejer.'''
    n = int((dn - d0 + timedelta(hours=30)).total_seconds() // (step_min * 60)) + 1
    t0 = d0 - timedelta(hours=30)
    ae = np.full(n, ae) if np.isscalar(ae) else np.asarray(ae, float)
    return np.array([[t0 + timedelta(minutes=step_min * k), float(ae[k]), f107] for k in range(n)], dtype=object)


def test_scherliess_fejer_shapes_and_sums():
    from mpyricalspace import models
    d0, dn = datetime(2015, 3, 17), datetime(2015, 3, 18)
    idx = _sf_indices(d0, dn)
    tt = np.array([d0 + timedelta(minutes=30 * k) for k in range(48)], dtype='datetime64[s]')
    ds = models.scherliess_fejer(tt, lon=-77., indices=idx)
    assert set(ds.data_vars) == {"qvdrift", "prompt", "dynamo", "dvdrift", "vdrift"}
    assert ds.qvdrift.size == 48
    assert np.array_equal(ds.time.values, tt)                       # evaluated at the given points
    assert np.allclose(ds.dvdrift.data, ds.prompt.data + ds.dynamo.data, atol=1e-6)
    assert np.allclose(ds.vdrift.data, ds.qvdrift.data + ds.dvdrift.data, atol=1e-6)
    assert np.all(np.abs(ds.qvdrift.data) < 100)                    # JRO quiet drift, m/s


def test_scherliess_fejer_prompt_responds_to_ae_step():
    from mpyricalspace import models
    d0, dn = datetime(2015, 3, 17, 6), datetime(2015, 3, 17, 18)
    n = int((dn - d0 + timedelta(hours=30)).total_seconds() // 900) + 1
    onset = int((30 + 3) * 4)                                       # AE jumps at ~d0 + 3 h
    ae = np.full(n, 80.); ae[onset:] = 900.
    idx = _sf_indices(d0, dn, ae=ae)
    tt = np.array([d0 + timedelta(minutes=15 * k) for k in range(49)], dtype='datetime64[s]')
    quiet_ae = models.scherliess_fejer(tt, lon=-77., indices=_sf_indices(d0, dn, ae=80.))
    storm    = models.scherliess_fejer(tt, lon=-77., indices=idx)
    assert np.allclose(quiet_ae.prompt.data, 0.0, atol=1e-9)        # constant AE -> no prompt penetration
    assert np.abs(storm.prompt.data).max() > 1.0                    # an AE jump drives it


def test_scherliess_fejer_irregular_times_and_deprecation():
    from mpyricalspace import models
    d0 = datetime(2015, 3, 17, 12)
    idx = _sf_indices(d0, d0 + timedelta(hours=1))
    odd = np.array([d0, d0 + timedelta(minutes=7), d0 + timedelta(minutes=41)], dtype='datetime64[s]')
    with pytest.warns(DeprecationWarning):
        ds = models.scherliess_fejer(odd, freq=timedelta(minutes=15), lon=-77.,
                                     indices=idx, half_resolution=True)
    assert np.array_equal(ds.time.values, odd)


def test_ppeefm1_rteef_local():
    '''the local RTEEF model runs from the bundled ACE data and returns 5-min rows.'''
    from mpyricalspace import models
    d0, dn = datetime(2004, 11, 9, 0), datetime(2004, 11, 10, 0)
    data = models._rteef(d0, dn, lon=-77.)
    assert data.shape == (289, 4)                              # 24 h at 5 min, + endpoint
    qef, tef, ppef = data[:, 1].astype(float), data[:, 2].astype(float), data[:, 3].astype(float)
    assert np.allclose(tef, qef + ppef)
    assert np.nanmax(np.abs(ppef)) < 20.0                      # prompt penetration, mV/m
    assert 0.0 <= np.nanmax(qef) < 2.0                         # daytime climatology, mV/m


def test_ppeefm1_out_of_range():
    from mpyricalspace import models
    with pytest.raises(ValueError):
        models._rteef(datetime(2015, 1, 1), datetime(2015, 1, 2))


def test_ppeefm1_nativelypackaged_indices():
    '''nativelypackaged_indices=True runs off the bundled SPIDR F107.txt instead of
       DataManager; ppef (independent of F10.7) is unaffected either way.'''
    from mpyricalspace import models
    d0, dn = datetime(2004, 11, 9, 0), datetime(2004, 11, 10, 0)
    default = models._rteef(d0, dn, lon=-77.)
    native  = models._rteef(d0, dn, lon=-77., nativelypackaged_indices=True)
    assert native.shape == default.shape == (289, 4)
    assert np.allclose(native[:, 3].astype(float), default[:, 3].astype(float))   # ppef


def test_msis_nativelypackaged_indices():
    '''nativelypackaged_indices=True pulls F10.7/ap from pymsis's own cache (CelesTrak)
       instead of DataManager; both paths return finite, physically sane density.'''
    o = Empirical()
    o.set_time([datetime(2013, 6, 15, 12)])
    o.set_location(-12., -77., 400.)
    default = o.run_msis()
    native  = o.run_msis(nativelypackaged_indices=True)
    for ds in (default, native):
        assert np.isfinite(float(ds.rho.data)) and float(ds.rho.data) > 0


def test_manoj_maus_ppeefm1(monkeypatch):
    from mpyricalspace import models
    d0, dn = datetime(2004, 11, 9, 12), datetime(2004, 11, 9, 18)
    tt = np.array([d0 + timedelta(hours=i) for i in range(7)], dtype='datetime64[s]')
    monkeypatch.setattr(models, 'scherliess_fejer', lambda *a, **k: xr.Dataset(
        {'qvdrift': ('time', np.full(7, 20.))},
        coords={'time': tt.astype('datetime64[s]').astype(datetime)}))
    ds = models.manoj_maus(tt, d0, dn, timedelta(hours=1), lon=-77., nativelypackaged_code=True)
    assert set(ds.data_vars) == {"mV_to_ms", "vdrift_manoj", "prompt_manoj", "qvdrift"}
    assert np.all(np.isfinite(ds.vdrift_manoj.data))


def test_run_methods_are_pure(obj):
    obj.run_igrf()
    obj.run_rocsat(f107s=np.full(obj.time.size, 150.))
    assert not hasattr(obj, "data")                                       # no state accumulation


def test_ndgrid_requires_time():
    from mpyricalspace._grid import flatten_ndgrid
    with pytest.raises(ValueError):
        flatten_ndgrid(lon=[0.0], add_slt=True)


def test_ndgrid_aligned_vs_meshgrid():
    from mpyricalspace._grid import flatten_ndgrid
    _, _, shp_aligned = flatten_ndgrid(lat=[0., 1., 2.], lon=[10., 11., 12.])   # same length -> hypercube
    _, _, shp_mesh = flatten_ndgrid(lat=[0., 1., 2.], lon=[10., 11.])           # different -> meshgrid
    assert shp_aligned == (3,)
    assert shp_mesh == (3, 2)
