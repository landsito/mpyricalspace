'''Backend selection + the duckdb query path (no web, no mongo).'''
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import pytest
from mpyricalspace.DataManager import DataManager


@pytest.fixture
def dm(tmp_path):
    d = DataManager(local=True, path=str(tmp_path), auto_fetch=False)
    d._write('daily', pd.DataFrame({'datetime': pd.date_range('2024-05-01', '2024-05-31'),
                                    'f107': np.linspace(200, 210, 31), 'f107a': 180.}))
    d._write('h3', pd.DataFrame({'datetime': pd.date_range('2024-05-01', '2024-05-31 21:00', freq='3h'),
                                 'ap': 100., 'kp': 5.0}))
    d._views()
    return d


def test_forces_local():
    assert DataManager(local=True, auto_fetch=False).db is None


def test_get_floors_to_cadence(dm):
    # 20:07 on the 11th -> f107 looked up at day floor, ap at 3h floor
    r = dm.get([datetime(2024, 5, 11, 20, 7)], ['f107', 'ap', 'kp'])
    assert r.shape == (1, 4)
    assert abs(float(r[0, 1]) - np.linspace(200, 210, 31)[10]) < 1e-6      # f107 of the 11th
    assert float(r[0, 2]) == 100.0 and float(r[0, 3]) == 5.0


def test_get_values_and_range(dm):
    v = dm.get_values([datetime(2024, 5, 3), datetime(2024, 5, 4)], 'f107')
    assert v.shape == (2,) and np.all(np.isfinite(v))
    rng = dm.get_range(datetime(2024, 5, 3, 0), datetime(2024, 5, 3, 2), ['ap'])
    assert rng.shape == (121, 2)                                          # minute grid, inclusive


def test_get_history(dm):
    h = dm.get_history([datetime(2024, 5, 11, 12)], 'ap', [-4, 0], timedelta(hours=3))
    assert h.shape == (1, 5) and np.all(h == 100.0)


def test_missing_param_is_nan(dm):
    r = dm.get([datetime(2024, 5, 11)], ['ae_15min'])                     # never written
    assert np.isnan(float(r[0, 1]))


def test_management_blocked_without_writes():
    dm = DataManager(uri="mongodb://127.0.0.1:1/", writes=False)          # not connected, never queried
    with pytest.raises(RuntimeError):
        dm.fetch('2024-01-01', '2024-01-02')


# ---------------------------------------------------------------- solar wind (OMNIWeb 1-min / 5-min)
SW_TEXT = """Selected parameters:
 1 Field magnitude average, nT
YYYY DOY HR MN      1       2       3       4       5      6        7     8      9
2024 132  2  3 9999.99 9999.99 9999.99 9999.99 99999.9 999.99 9999999. 99.99 999.99
2024 132  2  4   15.71   -8.10  -31.05    5.20   512.3   4.51  1234567.  3.10   -4.20
"""


@pytest.fixture
def sw(tmp_path):
    d = DataManager(local=True, path=str(tmp_path), auto_fetch=False)
    t1 = pd.date_range('2024-05-11 02:00', periods=3000, freq='min')
    d._write('sw1', pd.DataFrame({'datetime': t1, 'by': np.arange(3000.), 'bz': -5., 'vsw': 450.}))
    t5 = pd.date_range('2024-05-11 02:00', periods=100, freq='5min')
    d._write('sw5', pd.DataFrame({'datetime': t5, 'by_5min': np.arange(100.) * 10, 'bz_5min': np.nan}))
    d._views()
    return d


def test_sw_names_and_mongo_keys():
    import importlib
    M = importlib.import_module('mpyricalspace.DataManager')                       # the module, not the class of the same name
    assert M._KEY['by'] == M._KEY['by_5min'] == 'BY,_nT_(GSM)_(NOAA-OMNI)'          # the long OMNI field names used in the mongo collections
    assert M._CAD['bz'] == 'sw1' and M._CAD['bz_5min'] == 'sw5'


def test_sw_floors_to_native_cadence(sw):
    # 02:07 -> 1-min series at 02:07, 5-min series at the 02:05 sample
    r = sw.get([datetime(2024, 5, 11, 2, 7)], ['by', 'by_5min'])
    assert float(r[0, 1]) == 7.0 and float(r[0, 2]) == 10.0


def test_sw_missing_is_nan_not_none(sw):
    r = sw.get([datetime(2024, 5, 11, 2, 5), datetime(2024, 1, 1)], ['bz_5min', 'vsw'])
    assert all(isinstance(v, float) and np.isnan(v) for v in r[:, 1]) and np.isnan(float(r[1, 2]))
    assert float(r[0, 2]) == 450.0


def test_sw_history(sw):
    h = sw.get_history([datetime(2024, 5, 11, 2, 10)], 'by', [-3, 0], timedelta(minutes=1))
    assert h.tolist() == [[7.0, 8.0, 9.0, 10.0]]


def test_sw_dense_range_scan_matches_point_lookup(sw):
    # > _DENSE timestamps takes the BETWEEN path; result must equal the IN-list path
    dense = sw.get_range(datetime(2024, 5, 11, 2, 0), datetime(2024, 5, 13, 1, 59), ['by'])[:, 1].astype(float)
    assert dense.shape == (2880,) and dense[0] == 0.0 and dense[-1] == 2879.0
    few = sw.get([datetime(2024, 5, 11, 2, 0), datetime(2024, 5, 12, 12, 34)], ['by'])[:, 1].astype(float)
    assert few.tolist() == [dense[0], 24 * 60 + 10 * 60 + 34.]                              # 12 May 12:34 = 2074 min after 11 May 02:00


def test_parse_omni_fills_and_suffix():
    from mpyricalspace.DataManager import _parse_omni
    df = _parse_omni(SW_TEXT, 'min')
    assert list(df.columns) == ['datetime', 'bt', 'bx', 'by', 'bz', 'vsw', 'nsw', 'tsw', 'psw', 'esw']
    assert df.datetime.tolist() == [datetime(2024, 5, 11, 2, 3), datetime(2024, 5, 11, 2, 4)]     # DOY 132 of a leap year
    assert df.iloc[0, 1:].isna().all()                                                             # all fill codes -> NaN
    assert df.iloc[1].to_dict() | {'datetime': 0} == {'datetime': 0, 'bt': 15.71, 'bx': -8.10, 'by': -31.05, 'bz': 5.20,
                                                       'vsw': 512.3, 'nsw': 4.51, 'tsw': 1234567., 'psw': 3.10, 'esw': -4.20}
    assert _parse_omni(SW_TEXT, '5min').columns[3] == 'by_5min'


class _Col:
    'minimal pymongo collection stand-in'
    def __init__(self, docs=()):
        self.docs, self.ops, self.queries = list(docs), [], []
    def find(self, q, proj=None):
        self.queries.append((q, proj)); return iter(self.docs)
    def bulk_write(self, ops, ordered=False):
        self.ops += ops


def test_sw_mongo_routing_uses_long_keys():
    dm = DataManager(uri="mongodb://127.0.0.1:1/", writes=True)             # lazy client: never connects
    t = datetime(2024, 5, 11, 2, 5)
    dm.db = _Col()
    dm.sw = {'sw1': _Col(), 'sw5': _Col([{'datetime': t, 'BY,_nT_(GSM)_(NOAA-OMNI)': -31.0,
                                          'Speed,_km/s_(NOAA-OMNI)': 500.0}])}
    # read: 5-min collection, long keys projected, short names returned
    r = dm.get([datetime(2024, 5, 11, 2, 7)], ['by_5min', 'vsw_5min', 'bz_5min'])
    assert float(r[0, 1]) == -31.0 and float(r[0, 2]) == 500.0 and np.isnan(float(r[0, 3]))
    q, proj = dm.sw['sw5'].queries[0]
    assert q == {'datetime': {'$in': [t]}} and 'BY,_nT_(GSM)_(NOAA-OMNI)' in proj and not dm.db.queries
    # write: sw5 -> solar_wind.resolution_5min under the long key, one doc per sample, nothing forward-filled
    dm._write_mongo('sw5', pd.DataFrame({'datetime': [t], 'by_5min': [-31.0], 'bz_5min': [np.nan]}))
    assert len(dm.sw['sw5'].ops) == 1 and not dm.db.ops and not dm.sw['sw1'].ops
    op = dm.sw['sw5'].ops[0]
    assert op._filter == {'datetime': t} and op._doc == {'$set': {'BY,_nT_(GSM)_(NOAA-OMNI)': -31.0}}


def test_omni_text_clamps_to_available_range(monkeypatch):
    import requests
    from mpyricalspace.DataManager import _omni_text
    calls = []

    class R:
        def __init__(self, text): self.text = text

    def fake(url, params=None, timeout=None):
        calls.append(dict(params))
        return R("Error INVALID STOP DATE, correct range: 19810101 - 20260902" if params['end_date'] > '20260902' else SW_TEXT)
    monkeypatch.setattr(requests, 'get', fake)
    text, upto = _omni_text(datetime(2026, 8, 1), datetime(2026, 9, 1) + timedelta(days=30), 'min')
    assert text == SW_TEXT and upto == datetime(2026, 9, 3) and calls[1]['end_date'] == '20260902'

    monkeypatch.setattr(requests, 'get', lambda *a, **k: R("Error INVALID START DATE, correct range: 19810101 - 20260902"))
    text, upto = _omni_text(datetime(2026, 9, 15), datetime(2026, 10, 1), 'min')
    assert text is None and upto == datetime(2026, 9, 3)


def test_sw_fetch_on_mongo_backend_writes_solar_wind_only(monkeypatch):
    import importlib
    M = importlib.import_module('mpyricalspace.DataManager')
    monkeypatch.setattr(M, '_omni_text', lambda d0, dn, res: (SW_TEXT, None))
    dm = DataManager(uri="mongodb://127.0.0.1:1/", writes=True)
    dm.db, dm.sw = _Col(), {'sw1': _Col(), 'sw5': _Col()}
    dm.fetch('2024-05-01', '2024-05-02', include=['sw1', 'sw5'])            # must not touch duckdb views on mongo
    assert not dm.db.ops and len(dm.sw['sw1'].ops) == 1 and len(dm.sw['sw5'].ops) == 1
    assert dm.sw['sw1'].ops[0]._filter == {'datetime': datetime(2024, 5, 11, 2, 4)}      # native stamp, not forward-filled
