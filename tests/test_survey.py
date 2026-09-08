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
    assert {"igrf", "msis", "hwm", "iri", "hltwim", "sf", "rocsat", "eej", "eef", "manoj"} == set(d)
    assert d["iri"]["kind"] == "global" and d["eej"]["kind"] == "equator"
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
