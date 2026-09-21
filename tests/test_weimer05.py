'''Weimer 2005: compiled driver, dipole tilt and the python wrapper (no network, no mongo).'''
from datetime import datetime
import numpy as np
import pandas as pd
import pytest
from mpyricalspace import models
from mpyricalspace.DataManager import DataManager

# Reference values: the reference implementation's own test run (by=0, bz=-5, tilt=0, 450 km/s, 9 cm^-3), which prints 4 digits.
# (mlat, mlt, epot [kV], fac [uA/m2])
REF = [(90.0, 0., -14.35, 0.2876), (84.3344382773, 6., 4.384, -0.2093), (80.4295344893, 12., -14.61, 0.08036),
       (74.1876988925, 18., -52.95, -0.4077), (69.6681589023, 3., 35.78, 0.4099), (64.7975706790, 9., 1.292, -0.1744)]
T0 = datetime(2024, 5, 11, 12)
DRV = dict(by=0., bz=-5., vsw=450., nsw=9., tilt=0.)


def test_matches_original_program():
    lat, mlt, e, f = map(np.array, zip(*REF))
    ds = models.weimer05([T0] * len(REF), lat, mlt, **DRV)              # aligned samples
    assert np.allclose(ds.epot.values, e, rtol=2e-3) and np.allclose(ds.fac.values, f, rtol=2e-3)


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
    assert np.array_equal(a.epot.values, b.epot.values) and float(a.vsw) == 500.
    five = models.weimer05(T0, [75., 80.], [6., 12.], res='5min', avg=None, tilt=0)  # the 5-min series
    assert (float(five.by), float(five.vsw)) == (2.0, 600.)


def test_default_drivers_are_the_20min_mean_of_the_5min_series(store):
    # the store's 5-min bz ramps so that a 20-min mean differs from the instantaneous value
    t5 = pd.date_range('2024-05-11 11:00', '2024-05-11 12:30', freq='5min')
    store._write('sw5', pd.DataFrame({'datetime': t5, 'by_5min': 2.0, 'bz_5min': np.arange(len(t5)) * 1.0,
                                      'vsw_5min': 600., 'nsw_5min': 5.}))
    store._views()
    ds = models.weimer05(T0, [80.], [12.], tilt=0)                                # no res / avg: the defaults
    i = list(t5).index(pd.Timestamp(T0))
    assert float(ds.bz) == pytest.approx(np.mean([i - 3, i - 2, i - 1, i]))         # samples 11:45, 11:50, 11:55, 12:00
    explicit = models.weimer05(T0, [80.], [12.], res='5min', avg=20, tilt=0)
    assert float(explicit.bz) == float(ds.bz) and float(models.weimer05(T0, [80.], [12.], avg=None, tilt=0).bz) == float(i)


def test_trailing_average(store):
    ds = models.weimer05(T0, [80.], [12.], res='1min', avg=10, tilt=0)
    expect = np.mean(store.get_history([T0], 'bz', [-9, 0], pd.Timedelta(minutes=1).to_pytimedelta())[0])
    assert float(ds.bz) == pytest.approx(expect)


def test_lag_shifts_the_averaging_window(store):
    ds = models.weimer05(T0, [80.], [12.], res='1min', avg=5, lag=10, tilt=0)  # mean of the 5 min ending 10 min before
    expect = np.mean(store.get_history([T0], 'bz', [-14, -10], pd.Timedelta(minutes=1).to_pytimedelta())[0])
    assert float(ds.bz) == pytest.approx(expect)
    one = models.weimer05(T0, [80.], [12.], res='1min', avg=None, lag=10, tilt=0)   # no avg: the value 10 min earlier
    assert float(one.bz) == pytest.approx(float(store.get_values([T0 - pd.Timedelta(minutes=10).to_pytimedelta()], 'bz')[0]))


def test_no_data_gives_nan_not_error(store):
    ds = models.weimer05(datetime(2024, 5, 11, 3), [80.], [12.], tilt=0)          # nothing stored then
    assert np.isnan(ds.epot.values).all() and np.isnan(float(ds.by))


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
