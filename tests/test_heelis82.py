"""Heelis convection potential as NCAR TIE-GCM 2.0 computes it: compiled driver (both variants), the Kp -> cross-cap
potential relation, and the python wrapper (no network, no mongo)."""
from datetime import datetime
import numpy as np
import pandas as pd
import pytest
from mpyricalspace import models
from mpyricalspace.DataManager import DataManager

# Reference values: the UNTOUCHED src/heelis.F of TIE-GCM 2.0, driven with the parameters of its sub aurora_cons
# transcribed independently (cross-cap potential cp [kV], By [nT], latitude, MLT -> potential [V]).
# The unmodified file was run over 3600 points (3 cp x 6 By x both hemispheres x 10 lat x 10 MLT); every 89th is kept.
REF = [
       (15, -15, 30.5, 0, -0.004156691102569212),
       (15, -15, 89.9, 23.5, 3511.373749012095),
       (15, -15, -85, 21, -5555.327005616947),
       (15, -5, 80, 18, -8308.412436593691),
       (15, -5, -77, 15, -3574.994951999995),
       (15, 0, 76, 12, -1933.195196736448),
       (15, 0, -70, 11, -540.6231316253234),
       (15, 3, 60, 9, -2.215385629618029),
       (15, 3, -45, 6, 0.2474276355020202),
       (15, 10, 30.5, 3, 0.01259074699686899),
       (15, 10, 90, 0, -6783.305325958061),
       (15, 10, -85, 23.5, 330.4994905093058),
       (15, 7, 80, 21, -3904.242652736164),
       (15, 7, -77, 18, -5566.160987106784),
       (45, -15, 76, 15, -15905.10863903807),
       (45, -15, -70, 12, -7933.871310286884),
       (45, -5, 60, 11, -75.57840750593367),
       (45, -5, -45, 9, -1.618490005832312),
       (45, 0, 30.5, 6, 0.3736835825548656),
       (45, 0, 90, 3, -8887.642825741616),
       (45, 0, -89.9, 0, -8787.624694850942),
       (45, 3, 80, 23.5, -447.9051644453994),
       (45, 3, -77, 21, -17726.59285206644),
       (45, 10, 76, 18, -17615.77785548854),
       (45, 10, -70, 15, -7765.499329179003),
       (45, 7, 60, 12, -556.4803341389571),
       (45, 7, -45, 11, -0.216714827166584),
       (134, -15, 30.5, 9, 8.056068689976238),
       (134, -15, 90, 6, 26283.99867266642),
       (134, -15, -89.9, 3, -71241.8506612687),
       (134, -5, 85, 0, 2074.579326041539),
       (134, -5, -77, 23.5, 3170.406607518198),
       (134, 0, 76, 21, -46105.2185534286),
       (134, 0, -70, 18, -52370.05431230779),
       (134, 3, 60, 15, -10676.67463028588),
       (134, 3, -45, 12, -110.276684025877),
       (134, 10, 30.5, 11, -11.78411902137034),
       (134, 10, 90, 9, -55979.08654795062),
       (134, 10, -89.9, 6, 8398.15854016615),
       (134, 7, 85, 3, -24396.43122587773),
       (134, 7, -80, 0, 2342.835441643754),
]
T0 = datetime(2024, 5, 11, 12)
KEYS = ('theta0', 'psim', 'psie', 'pcen', 'phid', 'phin', 'phidp', 'phidm', 'phinp', 'phinm', 'offc', 'dskofc',
        'r1', 'r2', 'thetac', 'dtheta1', 'dtheta2')


def h(mlat, mlt, **kw):
    kw.setdefault('by', 0.)
    return models.heelis82(T0, mlat, mlt, **kw)


# ---------------------------------------------------------------------------- the TIE-GCM reference
def test_tiegcm_variant_matches_the_untouched_tiegcm_file():
    cp, by, lat, mlt, v = map(np.array, zip(*REF))
    got = np.array([float(h([la], [m], cp=c, by=b).epot) for c, b, la, m in zip(cp, by, lat, mlt)]) * 1e3
    assert np.allclose(got, v, rtol=1e-9, atol=1e-6)


def test_paper_variant_is_not_the_default():
    a, b = h(78., 18., cp=45.), h(78., 18., cp=45., variant='paper')
    assert float(a.epot) == float(h(78., 18., cp=45., variant='tiegcm').epot) and float(a.epot) != float(b.epot)


# ------------------------------------------------------------------------------ Kp -> cross-cap potential
@pytest.mark.parametrize('kp', [0., 1., 3., 4.5, 6., 9.])
def test_kp_to_cp_is_the_tiegcm_relation(kp):
    ds = h(78., 18., kp=kp)
    assert float(ds.cp) == pytest.approx(15. + 15. * kp + 0.8 * kp ** 2)                # util.F ctpoten_from_kp
    assert float(ds.epot) == float(h(78., 18., cp=float(ds.cp)).epot)


def test_cp_wins_over_kp():
    assert float(h(78., 18., kp=9., cp=45.).cp) == 45.
    assert float(h(78., 18., kp=9., cp=45.).epot) == float(h(78., 18., cp=45.).epot)


def test_bad_kp_or_cp_gives_nan_and_a_warning():
    for kw in ({'kp': -0.5}, {'kp': 9.5}, {'cp': -3.}, {'cp': 0.}):
        with pytest.warns(RuntimeWarning):
            assert np.isnan(float(h(78., 18., **kw).epot))
    assert np.isnan(float(h(78., 18., kp=np.nan).epot))                                 # missing: NaN, no warning


# ----------------------------------------------------------------------------------- the pattern
def test_extremes_are_the_boundary_potentials_and_their_difference_is_the_cross_cap_potential():
    mlat, mlt = np.arange(40., 90.01, 0.1), np.arange(0., 24., 0.05)
    e = h(mlat, mlt, cp=45.).epot.values
    assert e.max() == pytest.approx(0.44 * 45., abs=0.01) and e.min() == pytest.approx(-0.56 * 45., abs=0.01)
    assert e.max() - e.min() == pytest.approx(45., abs=0.02)


def test_morning_cell_positive_evening_cell_negative():
    lat = 77.
    assert float(h(lat, 5., cp=45.).epot) > 15. and float(h(lat, 17., cp=45.).epot) < -20.
    assert float(h(lat, 5., cp=90.).epot) > float(h(lat, 5., cp=45.).epot)                 # stronger drive


def test_zero_at_or_below_30_degrees_and_finite_at_the_poles():
    ds = h([-30., -20., 0., 20., 30., 30.001], [12.], cp=45.).epot.values.ravel()
    assert (ds[:5] == 0.).all() and ds[5] != 0.
    assert np.isfinite(h([90., -90.], [0., 12.], cp=45.).epot.values).all()


def test_by_limits_of_tiegcm():
    assert float(h(78., 18., cp=45., by=30.).epot) == float(h(78., 18., cp=45., by=7.).epot)
    assert float(h(78., 18., cp=45., by=-40.).epot) == float(h(78., 18., cp=45., by=-11.).epot)
    assert float(h(78., 18., cp=45., by=5.).epot) != float(h(78., 18., cp=45., by=0.).epot)


@pytest.mark.parametrize('variant', ['tiegcm', 'paper'])
def test_southern_hemisphere_is_the_northern_with_by_reversed(variant):
    for by in (4., -6.):
        s = h([-80., -70.], [3., 15.], cp=60., by=by, variant=variant)
        n = h([80., 70.], [3., 15.], cp=60., by=-by, variant=variant)
        assert np.allclose(s.epot.values, n.epot.values, rtol=1e-12, atol=1e-12)
    assert not np.allclose(h(-80., 3., cp=60., by=4.).epot, h(80., 3., cp=60., by=4.).epot)


def test_cube_shape_alignment_and_mixed_hemispheres():
    t = np.array(['2024-05-11T12:00', '2024-05-11T12:30'], dtype='datetime64[s]')
    lat, mlt = [80., -70., 75., -80.], [6., 18.]
    ds = models.heelis82(t, lat, mlt, cp=[45., 60.], by=[4., -2.])
    assert ds.epot.dims == ('time', 'mlat', 'mlt') and ds.epot.shape == (2, 4, 2)
    assert list(ds.cp.values) == [45., 60.] and list(ds.by.values) == [4., -2.]
    for ti in (0, 1):
        for i, la in enumerate(lat):
            one = models.heelis82(t[ti], [la], mlt, cp=[45., 60.][ti], by=[4., -2.][ti])
            assert np.array_equal(ds.epot.values[ti, i], one.epot.values.ravel())
    al = models.heelis82(t, [80., -70.], [6., 18.], cp=[45., 60.], by=[4., -2.])           # aligned samples
    assert al.epot.dims == ('time',) and al.epot.values[0] == ds.epot.values[0, 0, 0]


def test_input_validation():
    with pytest.raises(ValueError):
        h(78., 18., cp=45., variant='heelis')
    with pytest.raises(ValueError):
        h(78., 18., cp=45., res='10min')
    with pytest.raises(ValueError):
        h(78., 18., cp=45., params={'theta': 12.})
    t = np.array(['2024-05-11T12:00', '2024-05-11T12:01'], dtype='datetime64[s]')
    with pytest.raises(ValueError):
        models.heelis82(t, [80.], [6.], cp=[1, 2, 3], by=0.)


def test_params_override_the_constants():
    base = h(np.arange(60., 89., 2.), [3., 15.], cp=45.).epot.values
    same = h(np.arange(60., 89., 2.), [3., 15.], cp=45., params={'psim': 0.44 * 45., 'r1': -2.6}).epot.values
    assert np.allclose(base, same, rtol=1e-12)
    big = h(np.arange(60., 89., 2.), [3., 15.], cp=45., params={'psim': 40.}).epot.values
    assert not np.allclose(base, big) and big.max() > 30. > base.max()                     # psim = 40 kV raises the morning cell


# ----------------------------------------------------------- the variant of the paper (Heelis et al. 1982)
PAPER = dict(theta0=15., psim=20., psie=-35., pcen=0., phid=11., phin=24., phidp=22.5, phidm=22.5, phinp=55., phinm=55.,
             offc=0., dskofc=0., r1=-4., r2=2., thetac=14., dtheta1=1., dtheta2=1.)      # Fig. 4 + typical values, see README


def _G(colat, th0=15., dth1=1., dth2=1., r1=-4., r2=2., thc=14.):
    "latitude function of the paper (eqs. 2-3) with the constants in closed form"
    d, th0, th1, th2, thc = np.pi / 180, th0 * np.pi / 180, (th0 + dth1) * np.pi / 180, (th0 - dth2) * np.pi / 180, thc * np.pi / 180
    d1, d2 = th1 - th0, th2 - th0
    b1 = d1 * (d1 - np.tan(th1) / r1)
    a1 = np.sqrt(1 - d1 ** 2 / b1) / (np.sin(th1) / np.sin(th0)) ** r1
    c0 = (np.sin(thc) / np.sin(th0)) ** r2
    hh = (np.sin(th2 + thc) / np.sin(th0)) ** r2 - c0
    dh = r2 * (np.sin(th2 + thc) / np.sin(th0)) ** (r2 - 1) * np.cos(th2 + thc) / np.sin(th0)
    b2 = d2 * (d2 - hh / dh)
    a2 = np.sqrt(1 - d2 ** 2 / b2) / hh
    t = colat * d
    if t >= th1:
        return a1 * (np.sin(t) / np.sin(th0)) ** r1
    if t > th0:
        return np.sqrt(1 - (t - th0) ** 2 / b1)
    if t > th2:
        return np.sqrt(1 - (t - th0) ** 2 / b2)
    return a2 * ((np.sin(t + thc) / np.sin(th0)) ** r2 - c0)


def _F(mlt, colat, th0=15., pm=20e3, pe=-35e3):
    "local-time function of the paper (eqs. 4, 5a, 5b) for the Fig. 4 constants, in V; phi = hour angle clockwise from noon"
    ph = np.mod(-(mlt - 12.) * 15., 360.)
    fd, fn, x = 15., 180., ((colat - th0) / th0) ** 2
    ii, aa = 22.5 + x * (90 - 22.5), 55. + x * (90 - 55.)
    l1, l2 = min(90, (fn - fd) / 2), min(90, 180 - (fn - fd) / 2)
    ip, im, am, ap = min(ii, l1), min(ii, l2), min(aa, l1), min(aa, l2)
    a, b = (pm + pe) / 2, (pm - pe) / 2
    w = ph - 360 if ph > 180 else ph
    if fd - im <= w <= fd + ip:
        return a + b * np.cos(np.pi * (w - fd - ip) / (im + ip))
    if fn - am <= ph <= fn + ap:
        return a + b * np.cos(np.pi * (ph - fn + am) / (am + ap))
    return pm if fd + ip < ph < fn - am else pe


@pytest.mark.parametrize('colat', [0.5, 3., 8., 13.9, 14.5, 15.0001, 15.5, 16., 20., 24.])
def test_paper_variant_is_G_times_F_of_the_paper(colat):
    mlt = np.arange(0., 24., 0.25)
    e = h(90. - colat, mlt, cp=45., variant='paper', params=PAPER).epot.values * 1e3           # V
    ref = np.array([_G(colat) * _F(m, colat) for m in mlt])
    assert np.abs(e - ref).max() < 0.5                                                          # of a 55 kV swing


def test_paper_variant_smooths_the_flow_reversal_that_tiegcm_leaves_as_a_kink():
    lat = np.arange(73., 80.001, 0.05)
    kw = dict(cp=45., params={'offc': 0., 'dskofc': 0.})                                       # colat = 90 - lat
    jump = {}
    for v in ('tiegcm', 'paper'):
        g = np.diff(h(lat, [18.], variant=v, **kw).epot.values.ravel()) / 0.05                # kV/deg on the dusk meridian
        jump[v] = np.abs(np.diff(g)).max()
    assert jump['tiegcm'] > 2.5 and jump['paper'] < 0.6                                        # 3.5 vs 0.36 kV/deg


def test_paper_variant_keeps_the_cross_cap_potential_and_is_zero_at_the_pole():
    mlat, mlt = np.arange(40., 89.96, 0.1), np.arange(0., 24., 0.05)
    e = h(mlat, mlt, cp=45., variant='paper').epot.values
    assert e.max() - e.min() == pytest.approx(45., abs=0.05)
    assert abs(float(h(90., 12., cp=45., variant='paper', params={'offc': 1e-6}).epot)) < 0.2  # G(0) = 0 (tiegcm: ~ -9 kV)


def test_paper_variant_rejects_constants_without_a_solution():
    for bad in ({'dtheta2': 20.}, {'dtheta1': 0.}, {'r1': 1.}, {'r2': -1.}):
        with pytest.raises(ValueError):
            h(70., 12., cp=45., variant='paper', params=bad)
    assert np.isfinite(float(h(70., 12., cp=45., params={'dtheta2': 20.}).epot))              # the tiegcm variant has no such need


# ------------------------------------------------------------------------------ drivers from the DataManager
@pytest.fixture
def store(tmp_path, monkeypatch):
    dm = DataManager(local=True, path=str(tmp_path), auto_fetch=False)
    t5 = pd.date_range('2024-05-11 11:00', '2024-05-11 12:30', freq='5min')
    dm._write('sw5', pd.DataFrame({'datetime': t5, 'by_5min': np.arange(len(t5)) * 1.0}))
    t1 = pd.date_range('2024-05-11 11:00', '2024-05-11 12:30', freq='min')
    dm._write('sw1', pd.DataFrame({'datetime': t1, 'by': 2.0}))
    t3 = pd.date_range('2024-05-11 00:00', '2024-05-11 21:00', freq='3h')                     # Kp = 0 ... 7, one per interval
    dm._write('h3', pd.DataFrame({'datetime': t3, 'kp': np.arange(len(t3)) * 1.0}))
    dm._views()
    monkeypatch.setattr(models, 'manager', lambda: dm)
    return dm


def test_kp_is_interpolated_between_interval_centres_like_tiegcm(store):
    kp = lambda t, **kw: float(models._kp_from_store([t], **kw)[0])
    assert kp(datetime(2024, 5, 11, 13, 30)) == 4.0                        # the centre of the 12-15 interval (Kp = 4)
    assert kp(datetime(2024, 5, 11, 12, 0)) == 3.5                         # half way between the centres 10:30 and 13:30
    assert kp(datetime(2024, 5, 11, 15, 0)) == 4.5
    assert kp(datetime(2024, 5, 11, 12, 0), interp=False) == 4.0           # the value of the interval
    assert kp(datetime(2024, 5, 11, 14, 59), interp=False) == 4.0
    assert np.isnan(kp(datetime(2024, 5, 10, 22, 0)))                      # nothing stored the day before


def test_default_drivers_kp_interpolated_by_20min_mean(store):
    ds = models.heelis82(T0, [78.], [18.])
    i = list(pd.date_range('2024-05-11 11:00', '2024-05-11 12:30', freq='5min')).index(pd.Timestamp(T0))
    assert float(ds.kp) == 3.5 and float(ds.by) == pytest.approx(np.mean([i - 3, i - 2, i - 1, i]))
    assert float(ds.cp) == pytest.approx(15. + 15. * 3.5 + 0.8 * 3.5 ** 2)
    explicit = models.heelis82(T0, [78.], [18.], cp=float(ds.cp), by=float(ds.by))
    assert float(explicit.epot) == float(ds.epot)


def test_by_recipes(store):
    inst = models.heelis82(T0, [78.], [18.], res='1min', avg=None)
    assert float(inst.by) == 2.0
    five = models.heelis82(T0, [78.], [18.], res='5min', avg=None)
    i = list(pd.date_range('2024-05-11 11:00', '2024-05-11 12:30', freq='5min')).index(pd.Timestamp(T0))
    assert float(five.by) == float(i)
    assert float(models.heelis82(T0, [78.], [18.], avg=10, lag=5).by) == pytest.approx(i - 1.5)   # samples i-2, i-1


def test_cp_given_needs_no_kp_and_passed_values_are_not_averaged(store):
    ds = models.heelis82(datetime(2024, 5, 10, 3), [78.], [18.], cp=45., by=3.)         # no Kp stored for that day
    assert np.isfinite(float(ds.epot)) and float(ds.by) == 3.0


def test_no_data_gives_nan_not_error(store):
    ds = models.heelis82(datetime(2024, 5, 10, 3), [78.], [18.])
    assert np.isnan(float(ds.epot)) and np.isnan(float(ds.kp))


def test_facade_run_heelis82(store):
    from mpyricalspace.Predictor import Empirical
    obj = Empirical()
    obj.set_time([T0])
    ds = obj.run_heelis82(np.arange(50., 90.1, 5.), np.arange(0., 24., 6.))
    assert ds.epot.dims == ('mlat', 'mlt') and np.isfinite(ds.epot.values).all()
    assert (ds.epot.values == models.heelis82(T0, np.arange(50., 90.1, 5.), np.arange(0., 24., 6.)).epot.values).all()
