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
