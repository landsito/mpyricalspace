'''survey.run_track / run_grid -- multi-model evaluation over a track or a grid.'''
import numpy as np
import pytest
from datetime import datetime, timedelta

from mpyricalspace import survey


def _orbit(n=40, alt_km=400.0, year=2004):
    t0 = datetime(year, 3, 20, 12, 0)                                # noon -> daytime near the equator
    times = [t0 + timedelta(minutes=1.6 * i) for i in range(n)]
    lats = 25.0 * np.sin(np.linspace(0, 4 * np.pi, n))               # +-25 deg, crosses the equator
    lons = np.linspace(-120.0, 120.0, n)
    alts = np.full(n, alt_km)
    return times, lats, lons, alts


def test_model_domains():
    d = survey.model_domains()
    assert {"igrf", "msis", "hwm", "iri", "hltwim", "weimer05", "sf", "rocsat", "eej", "eef", "manoj"} == set(d)
    assert d["iri"]["kind"] == "global" and d["eej"]["kind"] == "equator"
    assert d["weimer05"]["kind"] == "high_lat" and d["weimer05"]["altitude_km"] == (0.0, 2000.0)
    assert d["sf"]["maglat_deg"] == 2.5 and d["eej"]["altitude_km"] == (90.0, 130.0)


def test_run_track_global_vs_equatorial():
    times, lats, lons, alts = _orbit()
    ds = survey.run_track(times, lats, lons, alts)

    assert ds.sizes["time"] == len(times)
    assert set(ds.attrs["models"]) >= {"igrf", "msis", "hwm", "iri", "sf", "eef", "manoj"}
    assert ds.attrs["skipped"] == {}

    # global: finite everywhere
    assert np.isfinite(ds.igrf_B.data).all()
    assert np.isfinite(ds.iri_NMF2.data).all()

    dip = ds.dip_lat.data
    # sf quiet drift: only within the +-2.5 deg cap, nowhere beyond
    assert np.isfinite(ds.sf_qvdrift.data[np.abs(dip) > 3.0]).sum() == 0
    assert np.isfinite(ds.sf_qvdrift.data[np.abs(dip) < 2.0]).any()
    # EEJ (E-region 90-130 km) is all-NaN on an F-region track
    assert not np.isfinite(ds.eej.data).any()


def test_run_track_eej_needs_e_region():
    times, lats, lons, _ = _orbit(alt_km=None)
    low = survey.run_track(times, lats, lons, alts=110.0, models=["eej", "sf"])
    dip = low.dip_lat.data
    assert np.isfinite(low.eej.data[np.abs(dip) < 3.0]).any()        # EEJ finite near the dip equator at 110 km
    assert not np.isfinite(low.sf_qvdrift.data).any()                # SF (F region) is NaN at 110 km


def test_run_track_manoj_date_gated():
    times, lats, lons, alts = _orbit(n=8, year=2015)
    ds = survey.run_track(times, lats, lons, alts, models=["sf", "manoj"])
    assert "sf" in ds.attrs["models"]
    assert "manoj" in ds.attrs["skipped"] and "2001-2007" in ds.attrs["skipped"]["manoj"]


def test_run_track_subset_and_altitude_default():
    times, lats, lons, _ = _orbit(12, alt_km=None)
    with pytest.warns(UserWarning):
        ds = survey.run_track(times, lats, lons, models=["igrf", "iri"])          # no alts -> default + warn
    assert ds.attrs["models"] == ["igrf", "iri"]
    assert not any(k.startswith("eej") for k in ds.data_vars)
    assert np.isfinite(ds.iri_NMF2.data).all()


def test_run_track_unknown_model_is_skipped():
    times, lats, lons, alts = _orbit(6)
    ds = survey.run_track(times, lats, lons, alts, models=["igrf", "not_a_model"])
    assert "not_a_model" in ds.attrs["skipped"] and ds.attrs["models"] == ["igrf"]


def test_run_grid_dims():
    ds = survey.run_grid(times=[datetime(2004, 6, 21, 12)],
                         lats=np.arange(-40.0, 41.0, 20.0), lons=np.arange(-135.0, 136.0, 90.0),
                         alts=[350.0], models=["igrf", "sf"])
    assert ds.igrf_B.dims == ("lat", "lon") and ds.igrf_B.shape == (5, 4)
    assert np.isfinite(ds.igrf_B.data).all()
    fin = np.isfinite(ds.sf_qvdrift.data)
    assert fin[ds.lat.data == 0].any()
    assert not fin[np.abs(ds.lat.data) == 40].any()


# ------------------------------------------------------------------------------------------- weimer05
DRV = dict(by=0., bz=-5., vsw=450., nsw=9., tilt=0.)         # explicit drivers: no index store, no network


def _polar():
    t0 = datetime(2024, 5, 12, 4, 0)
    times = [t0 + timedelta(minutes=k) for k in range(6)]
    lats = [65., 70., 80., -70., 20., 72.]                   # 3 north, 1 south, 1 too low, 1 above 2000 km
    lons = [-100., -98., -96., 120., 10., -94.]
    alts = [400., 400., 400., 400., 400., 3000.]
    return times, lats, lons, alts


def test_weimer05_default_run_at_low_latitude_is_nan_and_needs_nothing():
    times, lats, lons, alts = _orbit(n=8)                    # +-25 deg: the auroral cap is out of reach
    ds = survey.run_track(times, lats, lons, alts, models=["weimer05"])
    assert ds.attrs["models"] == ["weimer05"] and ds.attrs["skipped"] == {}
    assert not np.isfinite(ds.weimer05_epot.data).any() and not np.isfinite(ds.weimer05_fac.data).any()


def test_run_track_weimer05_matches_a_direct_call():
    import aacgmv2                                            # a package dependency: missing = a broken install
    from mpyricalspace import models
    times, lats, lons, alts = _polar()
    ds = survey.run_track(times, lats, lons, alts, models=["weimer05"], weimer_kw=DRV)
    assert ds.attrs["models"] == ["weimer05"] and ds.attrs["skipped"] == {}
    assert {"weimer05_epot", "weimer05_fac", "weimer05_mlat", "weimer05_mlt", "weimer05_tilt"} <= set(ds.data_vars)

    # the same samples converted one by one (exact time) and passed to models.weimer05 as aligned samples
    ok = [0, 1, 2, 3]
    ml, mt = [], []
    for i in ok:
        a, b, _ = aacgmv2.convert_latlon(lats[i], lons[i], alts[i], times[i], method_code="G2A")
        ml.append(a)
        mt.append(float(np.asarray(aacgmv2.convert_mlt(b, times[i])).reshape(-1)[0]))
    direct = models.weimer05([times[i] for i in ok], ml, mt, **DRV)
    assert np.allclose(ds.weimer05_mlat.data[ok], ml, atol=1e-2) and np.allclose(ds.weimer05_mlt.data[ok], mt, atol=1e-2)
    assert np.allclose(ds.weimer05_epot.data[ok], direct.epot.values, atol=0.1, equal_nan=True)
    assert np.isfinite(ds.weimer05_epot.data[:3]).all() and np.isfinite(ds.weimer05_fac.data[:3]).all()
    assert np.isfinite(ds.weimer05_epot.data[3])                             # southern hemisphere: the mirrored model
    assert not np.isfinite(ds.weimer05_epot.data[4:]).any()                  # 20 deg latitude; and above 2000 km


def test_run_track_weimer05_single_sample():
    import aacgmv2  # noqa: F401
    ds = survey.run_track([datetime(2024, 5, 12, 4, 0)], [75.], [-100.], [400.], models=["weimer05"], weimer_kw=DRV)
    assert ds.sizes["time"] == 1 and np.isfinite(ds.weimer05_epot.data).all()


def test_run_grid_weimer05_agrees_with_run_track():
    import aacgmv2  # noqa: F401
    t = datetime(2024, 5, 12, 4, 0)
    lats, lons = np.arange(60.0, 90.1, 10.0), np.arange(-180.0, 180.0, 60.0)
    grid = survey.run_grid(times=[t], lats=lats, lons=lons, alts=[400.0], models=["weimer05"], weimer_kw=DRV)
    assert grid.weimer05_epot.dims == ("lat", "lon") and grid.weimer05_epot.shape == (4, 6)
    la, lo = np.meshgrid(lats, lons, indexing="ij")
    track = survey.run_track([t] * la.size, la.ravel(), lo.ravel(), [400.0] * la.size,
                             models=["weimer05"], weimer_kw=DRV)
    assert np.array_equal(grid.weimer05_epot.values.ravel(), track.weimer05_epot.values, equal_nan=True)
    assert np.isfinite(grid.weimer05_epot.values).sum() >= 15                # most of a 60-90 deg cap is inside


def test_weimer05_without_aacgmv2_is_skipped_with_a_message(monkeypatch):
    import sys
    monkeypatch.setitem(sys.modules, "aacgmv2", None)                        # `import aacgmv2` -> ImportError
    times, lats, lons, alts = _polar()
    ds = survey.run_track(times, lats, lons, alts, models=["weimer05"], weimer_kw=DRV)
    assert "aacgmv2" in ds.attrs["skipped"]["weimer05"] and ds.attrs["models"] == []

