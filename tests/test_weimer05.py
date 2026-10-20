'''Weimer 2005: compiled driver, dipole tilt and the python wrapper (no network, no mongo).'''
import os
import warnings
from datetime import datetime
import numpy as np
import pandas as pd
import pytest
from mpyricalspace import models
from mpyricalspace.DataManager import DataManager

T0 = datetime(2024, 5, 11, 12)
DRV = dict(by=0., bz=-5., vsw=450., nsw=9., tilt=0.)
REF_DIR = os.path.join(os.path.dirname(__file__), 'data', 'weimer05')


def _read_w05sc(name):
    'one output file of the original test_w05sc program: (mlt, mlat, values [mlat, mlt]); its fill (1e36) -> NaN'
    tok = open(os.path.join(REF_DIR, name)).read().split()
    nlon, nlat = int(tok[0]), int(tok[1])
    v = np.array(tok[2:], float)
    mlt, mlat, z = v[:nlon], v[nlon:nlon + nlat], v[nlon + nlat:].reshape(nlat, nlon)
    return mlt, mlat, np.where(z > 1e35, np.nan, z)


def test_matches_original_program():
    '''The reference is the 2005b (SCHA) model itself: the output of its original test program, test_w05sc.f90
       (`make test` in src/weimer05; by=0, bz=-5, tilt=0, 450 km/s, 9 cm^-3, 4 printed digits), on its full
       97 x 25 grid. Northern hemisphere only: the program does not mirror negative latitudes (the wrapper does,
       2005a sec. 5, tested below).'''
    mlt, mlat, e = _read_w05sc('wei05sc_epot_f90.dat')
    _, _, f = _read_w05sc('wei05sc_fac_f90.dat')
    north = mlat >= 0
    ds = models.weimer05(T0, mlat[north], mlt, **DRV)
    for got, ref, atol in ((ds.epot.values, e[north], 5e-3), (ds.fac.values, f[north], 5e-5)):
        assert np.array_equal(np.isnan(got), np.isnan(ref))                      # same auroral boundary
        assert np.isfinite(ref).sum() > 300
        ok = np.isfinite(ref)
        assert np.allclose(got[ok], ref[ok], rtol=6e-4, atol=atol)               # to the 4 printed digits


def test_mpot_laplacian_is_the_fac():
    '''2005b eq. (1): J_par = lap(psi)/mu0. The surface Laplacian of mpot (psi/mu0) on the 110 km sphere, by finite
       differences, must give `fac` with its sign and scale. Not point by point: the model computes fac from the
       series truncated at degree 10, whereas mpot keeps every degree.'''
    d = 0.25
    lat, mlt = np.arange(60., 86.01, d), np.arange(0., 24., d / 15.)       # the same angular step in both
    ds = models.weimer05(T0, lat, mlt, **DRV)
    psi, fac = ds.mpot.values * 1e3, ds.fac.values * 1e-6                   # A, A/m^2
    assert ds.mpot.attrs['units'] == 'kA' and np.isnan(psi[0]).any() and np.isfinite(psi).sum() > 1000
    r, h, cl = (6371.2 + 110.) * 1e3, np.radians(d), np.cos(np.radians(lat))[:, None]
    lap = (np.gradient(cl * np.gradient(psi, h, axis=0), h, axis=0) / cl
           + (np.roll(psi, -1, 1) - 2 * psi + np.roll(psi, 1, 1)) / h ** 2 / cl ** 2) / r ** 2
    ok = np.isfinite(lap) & (np.abs(fac) > 0.2e-6)
    ok[:4] = ok[-4:] = False                                                # one-sided differences at the grid edges
    assert ok.sum() > 1000
    assert np.corrcoef(lap[ok], fac[ok])[0, 1] > 0.9
    assert np.median(fac[ok] / lap[ok]) == pytest.approx(1., abs=0.1)


def test_nan_equatorward_of_the_boundary_raises_no_warning():
    'NaN equatorward of the auroral boundary is the normal output there, not something to warn about'
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ds = models.weimer05(T0, np.arange(-90., 90.1, 2.), np.arange(0., 24., 1.), **DRV)
    assert np.isnan(ds.epot.values).any() and np.isfinite(ds.epot.values).any()


def test_cube_shape_coords_and_attrs():
    t = np.array(['2024-05-11T12:00', '2024-05-11T12:01'], dtype='datetime64[s]')
    ds = models.weimer05(t, [70., 80., 85.], [0., 6., 12., 18.], by=[0, 1], bz=[-5, -6], vsw=450, nsw=9, tilt=0)
    assert ds.epot.dims == ('time', 'mlat', 'mlt') and ds.epot.shape == (2, 3, 4)
    assert ds.by.dims == ('time',) and list(ds.by.values) == [0., 1.] and (ds.vsw.values == 450.).all()
    assert ds.epot.attrs['units'] == 'kV' and ds.fac.attrs['units'] == 'uA/m^2'
    # each step is set up on its own drivers: step 1 must equal a stand-alone call with the same drivers
    one = models.weimer05(t[1], [70., 80., 85.], [0., 6., 12., 18.], by=1, bz=-6, vsw=450, nsw=9, tilt=0)
    assert np.array_equal(ds.epot.values[1], one.epot.values, equal_nan=True)


def test_nan_outside_cap_and_for_missing_driver():
    ds = models.weimer05([T0] * 3, [80., 75., 80.], [6., 12., 6.], by=[0, 0, np.nan], bz=-5, vsw=450, nsw=9, tilt=0)
    e, f = ds.epot.values, ds.fac.values
    assert np.isfinite(e[:2]).all() and np.isfinite(f[:2]).all()
    assert np.isnan(e[2]) and np.isnan(f[2])                                              # a NaN driver
    assert np.isnan(models.weimer05(T0, [50.], [6.], **DRV).epot.values).all()          # equatorward of the cap
    assert np.isnan(models.weimer05(T0, [-50.], [6.], **DRV).epot.values).all()


def test_southern_hemisphere_is_the_mirrored_north():
    # Weimer (2005a, sec. 5): the southern hemisphere uses the northern model with By and the tilt reversed
    d = dict(bz=-5., vsw=450., nsw=9.)
    s = models.weimer05(T0, [-80., -70.], [3., 15.], by=4., tilt=12., **d)
    n = models.weimer05(T0, [80., 70.], [3., 15.], by=-4., tilt=-12., **d)
    assert np.isfinite(s.epot.values).all()
    assert np.array_equal(s.epot.values, n.epot.values) and np.array_equal(s.fac.values, n.fac.values)
    # a cube with both hemispheres in any order equals the individual calls, and keeps the caller's mlat order
    lat = [80., -70., 75., -80.]
    cube = models.weimer05(np.array([T0, datetime(2024, 5, 11, 12, 30)], dtype='datetime64[s]'), lat, [6., 18.],
                           by=[4., -2.], tilt=[12., 5.], **d)
    for i, la in enumerate(lat):
        for ti in (0, 1):
            one = models.weimer05(cube.time.values[ti], [la], [6., 18.], by=[4., -2.][ti], tilt=[12., 5.][ti], **d)
            assert np.array_equal(cube.epot.values[ti, i], one.epot.values, equal_nan=True)


def test_input_validation():
    t = np.array(['2024-05-11T12:00', '2024-05-11T12:01'], dtype='datetime64[s]')
    with pytest.raises(ValueError):
        models.weimer05(t, [80.], [6.], by=[1, 2, 3], bz=-5, vsw=450, nsw=9, tilt=0)
    with pytest.raises(ValueError):
        models.weimer05(t, [80.], [6.], res='10min', **DRV)


def test_missing_data_files_raise(monkeypatch, tmp_path):
    monkeypatch.setattr(models, '_datadir', lambda name: str(tmp_path))
    with pytest.raises(FileNotFoundError):
        models.weimer05(T0, [80.], [6.], **DRV)


def test_extreme_drivers_warn_instead_of_crashing():
    with pytest.warns(RuntimeWarning):
        ds = models.weimer05(T0, [80.], [6.], by=0., bz=-60., vsw=1500., nsw=80., tilt=0.)
    assert np.isnan(ds.epot.values).all()


# ------------------------------------------------------------------ drivers from the DataManager
@pytest.fixture
def store(tmp_path, monkeypatch):
    dm = DataManager(local=True, path=str(tmp_path), auto_fetch=False)
    t = pd.date_range('2024-05-11 11:00', '2024-05-11 12:30', freq='min')
    dm._write('sw1', pd.DataFrame({'datetime': t, 'by': 1.0, 'bz': np.linspace(-8, -2, len(t)), 'vsw': 500., 'nsw': 7.}))
    t5 = pd.date_range('2024-05-11 11:00', '2024-05-11 12:30', freq='5min')
    dm._write('sw5', pd.DataFrame({'datetime': t5, 'by_5min': 2.0, 'bz_5min': -4., 'vsw_5min': 600., 'nsw_5min': 5.}))
    dm._views()
    monkeypatch.setattr(models, 'manager', lambda: dm)
    return dm


def test_drivers_from_datamanager(store):
    a = models.weimer05(T0, [75., 80.], [6., 12.], res='1min', avg=None, tilt=0)     # instantaneous 1-min values
    bz = float(store.get_values([T0], 'bz')[0])
    b = models.weimer05(T0, [75., 80.], [6., 12.], by=1., bz=bz, vsw=500., nsw=7., tilt=0)
    assert np.array_equal(a.epot.values, b.epot.values) and a.vsw.item() == 500.
    five = models.weimer05(T0, [75., 80.], [6., 12.], res='5min', avg=None, tilt=0)  # the 5-min series
    assert (five.by.item(), five.vsw.item()) == (2.0, 600.)


def test_default_drivers_are_the_20min_mean_of_the_5min_series(store):
    # the store's 5-min bz ramps so that a 20-min mean differs from the instantaneous value
    t5 = pd.date_range('2024-05-11 11:00', '2024-05-11 12:30', freq='5min')
    store._write('sw5', pd.DataFrame({'datetime': t5, 'by_5min': 2.0, 'bz_5min': np.arange(len(t5)) * 1.0,
                                      'vsw_5min': 600., 'nsw_5min': 5.}))
    store._views()
    ds = models.weimer05(T0, [80.], [12.], tilt=0)                                # no res / avg: the defaults
    i = list(t5).index(pd.Timestamp(T0))
    assert ds.bz.item() == pytest.approx(np.mean([i - 3, i - 2, i - 1, i]))         # samples 11:45, 11:50, 11:55, 12:00
    explicit = models.weimer05(T0, [80.], [12.], res='5min', avg=20, tilt=0)
    assert explicit.bz.item() == ds.bz.item() and float(models.weimer05(T0, [80.], [12.], avg=None, tilt=0).bz) == float(i)


def test_trailing_average(store):
    ds = models.weimer05(T0, [80.], [12.], res='1min', avg=10, tilt=0)
    expect = np.mean(store.get_history([T0], 'bz', [-9, 0], pd.Timedelta(minutes=1).to_pytimedelta())[0])
    assert ds.bz.item() == pytest.approx(expect)


def test_lag_shifts_the_averaging_window(store):
    ds = models.weimer05(T0, [80.], [12.], res='1min', avg=5, lag=10, tilt=0)  # mean of the 5 min ending 10 min before
    expect = np.mean(store.get_history([T0], 'bz', [-14, -10], pd.Timedelta(minutes=1).to_pytimedelta())[0])
    assert ds.bz.item() == pytest.approx(expect)
    one = models.weimer05(T0, [80.], [12.], res='1min', avg=None, lag=10, tilt=0)   # no avg: the value 10 min earlier
    assert one.bz.item() == pytest.approx(float(store.get_values([T0 - pd.Timedelta(minutes=10).to_pytimedelta()], 'bz')[0]))


def test_no_data_gives_nan_not_error(store):
    ds = models.weimer05(datetime(2024, 5, 11, 3), [80.], [12.], tilt=0)          # nothing stored then
    assert np.isnan(ds.epot.values).all() and np.isnan(ds.by.item())


# ------------------------------------------------------------------------------- dipole tilt
@pytest.mark.parametrize('when, expected', [('2024-06-21T18:00:00', 32.243), ('2000-12-21T06:00:00', -33.243),
                                            ('2024-03-20T12:00:00', 2.593), ('2031-09-01T00:00:00', 5.738)])
def test_dipole_tilt_against_astropy(when, expected):
    # expected: astropy Sun (ITRS) with the same IGRF-14 dipole; the wrapper's low-precision Sun is good to ~0.01 deg
    assert models.dipole_tilt(when)[0] == pytest.approx(expected, abs=0.02)


def test_dipole_tilt_seasons_and_shape():
    day = np.datetime64('2024-06-21') + np.arange(0, 1440, 10).astype('timedelta64[m]')
    june, dec = models.dipole_tilt(day), models.dipole_tilt(day + np.timedelta64(183, 'D'))
    assert june.shape == day.shape and 12 < june.min() and june.max() < 34 and dec.max() < -12   # summer / winter


def test_dipole_table_is_not_behind_the_newest_bundled_igrf():
    'dipole_tilt uses a copy of the IGRF dipole terms; this fails when a newer generation is vendored under src/igrf'
    import glob, os, re
    root = os.path.join(os.path.dirname(models.__file__), '..', 'src', 'igrf')
    dirs = sorted(glob.glob(os.path.join(root, 'igrf[0-9][0-9]')), key=lambda d: int(d[-2:]))
    if not dirs:
        pytest.skip('source tree not available (installed package)')
    epochs = {}
    for ln in open(glob.glob(os.path.join(dirs[-1], 'igrf*.f'))[0]):
        m = re.match(r'\s*data\s+(g[a-z]|g\d)\s*/(.*)', ln, re.I)
        if m:
            nums = re.findall(r'-?\d+\.?\d*', m.group(2))
            if int(float(nums[-1])) % 5 == 0:                       # an epoch row (the secular-variation row is not)
                epochs[int(float(nums[-1]))] = [float(x) for x in nums[:3]]
    for row in models._DIPOLE:
        assert int(row[0]) in epochs and np.allclose(row[1:], epochs[int(row[0])])
    assert max(epochs) <= models._DIPOLE[-1, 0], 'a newer IGRF epoch is bundled: refresh models._DIPOLE'


def test_igrf_latest_alias():
    a = models.igrf(T0, 10., 20., 300.)                            # default = 'latest'
    b = models.igrf(T0, 10., 20., 300., version=14)
    assert np.allclose(a.B.values, b.B.values) and np.allclose(models.igrf(T0, 10., 20., 300., version='latest').B.values, b.B.values)
    with pytest.raises(ValueError):
        models.igrf(T0, 10., 20., 300., version=15)
