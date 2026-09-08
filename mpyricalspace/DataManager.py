'''
Geophysical-index store for mpyricalspace.

Backend is chosen at runtime (nothing extra to install for either):
  * mongo   -> uri arg, $MPYRICALSPACE_MONGO_URI, config file
  * duckdb  -> local parquet under $MPYRICALSPACE_DATA_DIR (default ~/.cache/mpyricalspace)
`local=True` or $MPYRICALSPACE_LOCAL=1 forces duckdb.

Query API (identical output shapes on both backends, so Predictor.py does not care):
  dm = DataManager()
  dm.get(dates, ['f107','kp'])
  dm.get_range(d0, dn, ['ae_15min','f107'])
  dm.get_values(dates, 'f107')
  dm.get_history(dates, 'ap', [-20,0], timedelta(hours=3))

On the duckdb backend the store auto-fills any date span a model asks for (auto_fetch);
`python -m mpyricalspace {status,fetch,update,config}` drives it explicitly.

Download sources:
  kp/ap/f107/f107a  GFZ Potsdam        (source: potsdam|potsdam-legacy|noaa|canada|lisird)
  ae/ao/al/au       Kyoto WDC aeasy-cgi -> realtime-image scrape -> NOAA OMNI  (ae_source)
  dst               Kyoto WDC
  symh/asy          Kyoto WDC aeasy-cgi
'''
import os, glob, json, time, warnings, numpy as np
from datetime import datetime, timedelta
from pandas import to_datetime, date_range, DataFrame, Series, concat

DEFAULT_DIR = "~/.cache/mpyricalspace"
CONFIG      = os.path.expanduser("~/.config/mpyricalspace/config.json")
GFZ_START   = datetime(1932, 1, 1)
_env        = lambda k: os.environ.get(k) or None
_cfg        = lambda: json.load(open(CONFIG)) if os.path.exists(CONFIG) else {}
tofloat     = lambda w: np.nan if not w.strip() or '99999' in w else float(w)

# param -> native cadence ; param -> download group ; cadence -> left-aligned floor
_CAD = {'f107':'daily','f107adj':'daily','f107a':'daily','ap_daily':'daily','kp_daily':'daily','sn':'daily',
        'ap':'h3','kp':'h3',
        'ae_15min':'q15','ao_15min':'q15','au_15min':'q15','al_15min':'q15',
        'ae_hourly':'hourly','ao_hourly':'hourly','au_hourly':'hourly','al_hourly':'hourly','dst':'hourly',
        'ae':'min','ao':'min','au':'min','al':'min','symh':'min','symd':'min','asyh':'min','asyd':'min'}
_GRP = {**{p:'kpap' for p in ('f107','f107adj','f107a','ap_daily','kp_daily','sn','ap','kp')},
        **{p:'ae'   for p in ('ae','ao','al','au','ae_15min','ao_15min','au_15min','al_15min',
                              'ae_hourly','ao_hourly','au_hourly','al_hourly')},
        'dst':'dst', **{p:'asy' for p in ('symh','symd','asyh','asyd')}}
_FLR = {'daily' : lambda d: datetime(d.year, d.month, d.day),
        'h3'    : lambda d: d.replace(hour=d.hour // 3 * 3, minute=0, second=0, microsecond=0),
        'hourly': lambda d: d.replace(minute=0, second=0, microsecond=0),
        'q15'   : lambda d: d.replace(minute=d.minute // 15 * 15, second=0, microsecond=0),
        'min'   : lambda d: d.replace(second=0, microsecond=0)}
_dts = lambda x: [d.replace(second=0, microsecond=0) for d in to_datetime(np.atleast_1d(x)).to_pydatetime()]
_months = lambda t0, t1: [(datetime(y, m, 1), datetime(y + (m == 12), m % 12 + 1, 1))
                          for y in range(t0.year, t1.year + 1)
                          for m in range(1, 13) if datetime(y, m, 1) <= _FLR['daily'](t1)
                          and datetime(y + (m == 12), m % 12 + 1, 1) > _FLR['daily'](t0)]


class DataManager(object):

    def __init__(self, uri=None, path=None, local=False, auto_fetch=True, writes=None):
        cfg   = _cfg()
        local = local or _env("MPYRICALSPACE_LOCAL") not in (None, "0")
        uri   = None if local else (uri or _env("MPYRICALSPACE_MONGO_URI") or cfg.get("mongo_uri"))
        self.auto_fetch, self.db = auto_fetch, None
        self._cov, self._covp = {}, None
        # fetch/update may write to mongo only when explicitly allowed (protects a read-only production db)
        self.writes = writes if writes is not None else (_env("MPYRICALSPACE_MONGO_WRITE") not in (None, "0") or bool(cfg.get("mongo_writes")))
        if uri:
            import pymongo
            u = _env("MPYRICALSPACE_MONGO_USER") or cfg.get("mongo_user")
            p = _env("MPYRICALSPACE_MONGO_PASSWORD") or cfg.get("mongo_password")
            self.db = pymongo.MongoClient(uri, **({'username': u, 'password': p} if u else {})).geomagnetic_indices.resolution_1min
        if self.db is None:
            import duckdb
            self.path  = os.path.expanduser(path or _env("MPYRICALSPACE_DATA_DIR") or cfg.get("data_dir") or DEFAULT_DIR)
            self._covp = os.path.join(self.path, "_coverage.json")
            self._cov  = json.load(open(self._covp)) if os.path.exists(self._covp) else {}
            self.con   = duckdb.connect()
            self._views()

    def _views(self):
        self._have, self._cols = set(), {}
        for c in set(_CAD.values()):
            self.con.execute("DROP VIEW IF EXISTS %s" % c)
            g = os.path.join(self.path, c, '**', '*.parquet')
            if glob.glob(g):
                self.con.execute("CREATE VIEW %s AS SELECT * FROM read_parquet('%s', hive_partitioning=true, union_by_name=true)" % (c, g))
                self._have.add(c)
                self._cols[c] = {d[0] for d in self.con.execute("DESCRIBE %s" % c).fetchall()}

    # --------------------------------------------------------------------- query
    def _query(self, cad, params, keys):
        'floored datetime -> mapping {param: value}'
        if self.db is not None:
            return {d['datetime']: d for d in self.db.find({'datetime': {'$in': list(keys)}}, {p: 1 for p in ['datetime', *params]})}
        if cad not in getattr(self, '_have', ()) or not keys:
            return {}
        cols = [p for p in params if p in self._cols.get(cad, ())]
        if not cols:
            return {}
        q = "SELECT %s FROM %s WHERE datetime IN (%s)" % (",".join(['datetime', *cols]), cad, ",".join(['?'] * len(keys)))
        return {r[0]: dict(zip(['datetime', *cols], r)) for r in self.con.execute(q, list(keys)).fetchall()}

    def _lookup(self, params, dts):
        self._ensure(params, dts)
        out = [{} for _ in dts]
        for cad in {_CAD[p] for p in params}:
            ps, keys = [p for p in params if _CAD[p] == cad], [_FLR[cad](d) for d in dts]
            rows = self._query(cad, ps, sorted(set(keys)))
            for o, k in zip(out, keys):
                o.update({p: rows.get(k, {}).get(p, np.nan) for p in ps})
        return out

    def get(self, dates, params):
        dts = _dts(dates)
        return np.array([[d, *[o[p] for p in params]] for d, o in zip(dts, self._lookup(params, dts))], dtype=object)

    def get_values(self, dates, param):
        return self.get(dates, [param])[:, 1].astype(float)

    def get_range(self, d0, dn, params):
        return self.get(date_range(_dts(d0)[0], _dts(dn)[0], freq='min').to_pydatetime(), params)

    def get_history(self, dates, param, lags, res_lag):
        cad  = _CAD[param]
        dd   = _dts(dates)
        grid = [[_FLR[cad](d + res_lag * k) for k in range(lags[0], lags[1] + 1)] for d in dd]
        self._ensure([param], [m for row in grid for m in row])
        lut  = self._query(cad, [param], sorted({m for row in grid for m in row}))
        return np.array([[lut.get(m, {}).get(param, np.nan) for m in row] for row in grid], dtype=float)

    # ------------------------------------------------------------------ coverage
    def _covered(self, g, t0, t1):
        return any(a <= t0 and t1 <= b for a, b in self._cov.get(g, []))

    def _mark(self, g, t0, t1):
        iv = sorted(self._cov.get(g, []) + [[t0.isoformat(), t1.isoformat()]])
        m  = [iv[0]]
        for a, b in iv[1:]:
            if a <= m[-1][1]: m[-1][1] = max(m[-1][1], b)
            else: m.append([a, b])
        self._cov[g] = m
        if self._covp:
            os.makedirs(self.path, exist_ok=True)
            json.dump(self._cov, open(self._covp, 'w'), indent=0)

    def _ensure(self, params, dts):
        if self.db is not None or not self.auto_fetch or not dts:
            return
        t0, t1  = min(dts), max(dts)
        groups  = {_GRP[p] for p in params}
        touched = False
        fresh   = os.path.exists(self._covp) and (time.time() - os.path.getmtime(self._covp) < 21600)
        if 'kpap' in groups and not self._covered('kpap', t0.isoformat(), t1.isoformat()):
            touched = True
            try:
                a, b = t0 - timedelta(days=100), max(t1, datetime.utcnow()) + timedelta(days=1)
                self._fetch_kpap(a, b); self._mark('kpap', a, b)
            except Exception as e:
                warnings.warn("kp/ap auto-fetch failed (%s); check connectivity or fetch manually" % e)
        if 'ae' in groups:
            for m0, m1 in _months(t0, t1):
                stale = (datetime.utcnow() - m1).days < 60 and not fresh
                if stale or not self._covered('ae', m0.isoformat(), m1.isoformat()):
                    touched = True
                    try:
                        self._fetch_ae(m0, m1); self._mark('ae', m0, m1)
                    except Exception as e:
                        warnings.warn("AE auto-fetch for %s failed (%s)" % (m0.strftime('%Y-%m'), e))
        if touched:
            self._views()

    # ---------------------------------------------------------------- management
    def _writable(self):
        if self.db is not None and not self.writes:
            raise RuntimeError("mongo backend is read-only here; pass writes=True / $MPYRICALSPACE_MONGO_WRITE=1 "
                               "to let fetch/update populate it")

    def status(self):
        if self.db is not None:
            n = self.db.estimated_document_count()
            ext = self.db.find_one(sort=[('datetime', 1)]), self.db.find_one(sort=[('datetime', -1)])
            print("backend: mongo (%s docs, %s .. %s)%s" % (
                n, ext[0] and ext[0]['datetime'], ext[1] and ext[1]['datetime'],
                "" if self.writes else "   [read-only]"))
            return
        print("store: %s" % self.path)
        for c in sorted(set(_CAD.values())):
            names = ", ".join(p for p in _CAD if _CAD[p] == c)
            if c not in self._have:
                print("  %-7s empty                                    (%s)" % (c, names)); continue
            n, t0, t1 = self.con.execute("SELECT count(*), min(datetime), max(datetime) FROM %s" % c).fetchone()
            print("  %-7s %s .. %s  %9i rows  (%s)" % (c, t0, t1, n, names))
        for g, iv in self._cov.items():
            print("  coverage[%s]: %s" % (g, ", ".join("%s..%s" % (a[:10], b[:10]) for a, b in iv)))

    def fetch(self, d0, dn, kpap_source='potsdam', ae_source='auto', include=None):
        '''download every group (or `include` subset) over [d0, dn] and store it
           (parquet on the duckdb backend, forward-filled 1-min docs on mongo).
           include: any of  kpap ae dst asy'''
        self._writable()
        d0, dn = _dts(d0)[0], _dts(dn)[0]
        groups = include or ['kpap', 'ae', 'dst', 'asy']
        if 'kpap' in groups:
            self._fetch_kpap(d0, dn, kpap_source); self._mark('kpap', d0, dn)
        for m0, m1 in _months(d0, dn):
            if 'ae'  in groups: self._try(self._fetch_ae,  'ae',  m0, m1, ae_source)
            if 'dst' in groups: self._try(self._fetch_dst, 'dst', m0, m1)
            if 'asy' in groups: self._try(self._fetch_asy, 'asy', m0, m1)
        self._views()

    def update(self, since=None):
        'fill every group from the end of its coverage (or `since`) up to today'
        self._writable()
        dn = datetime.utcnow()
        if since is None and self.db is not None:
            last = self.db.find_one(sort=[('datetime', -1)])
            since = (last['datetime'] - timedelta(days=60)) if last else GFZ_START
        for g in ('kpap', 'ae', 'dst', 'asy'):
            iv = self._cov.get(g, [])
            d0 = since or (datetime.fromisoformat(iv[-1][1]) - timedelta(days=60) if iv else GFZ_START)
            try:    self.fetch(d0, dn, include=[g])
            except Exception as e: warnings.warn("update[%s] failed: %s" % (g, e))

    def _try(self, fn, g, m0, m1, *a):
        try:    fn(m0, m1, *a); self._mark(g, m0, m1)
        except Exception as e: warnings.warn("%s %s failed: %s" % (g, m0.strftime('%Y-%m'), e))

    _SPAN = {'daily': 1440, 'h3': 180, 'hourly': 60, 'q15': 15, 'min': 1}   # minutes a native-cadence row covers

    def _write(self, cad, df):
        'persist a native-cadence DataFrame (datetime + param columns) to the active backend'
        df = df.dropna(how='all', subset=[c for c in df.columns if c != 'datetime'])
        if df.empty:
            return
        if self.db is not None:
            return self._write_mongo(cad, df)
        for yr, chunk in df.groupby(to_datetime(df.datetime).dt.year):
            d = os.path.join(self.path, cad, "year=%i" % yr); os.makedirs(d, exist_ok=True)
            f = os.path.join(d, 'part.parquet')
            chunk = chunk.set_index('datetime')
            if os.path.exists(f):
                old   = self.con.execute("SELECT * FROM read_parquet('%s')" % f).df().set_index('datetime')
                chunk = chunk.combine_first(old)
            chunk = chunk.sort_index().reset_index()
            self.con.register('_w', chunk)
            self.con.execute("COPY _w TO '%s' (FORMAT parquet)" % f)
            self.con.unregister('_w')

    def _write_mongo(self, cad, df):
        'forward-fill each native-cadence row to 1-min docs and $set-upsert'
        from pandas import isna
        from pymongo import UpdateOne
        span, ops = self._SPAN[cad], []
        for _, r in df.iterrows():
            t0  = r['datetime'].to_pydatetime().replace(second=0, microsecond=0)
            doc = {k: (float(v) if hasattr(v, 'dtype') else v) for k, v in r.items() if k != 'datetime' and not isna(v)}
            if not doc:
                continue
            ops += [UpdateOne({'datetime': t0 + timedelta(minutes=i)}, {'$set': doc}, upsert=True) for i in range(span)]
            if len(ops) >= 5000:
                self.db.bulk_write(ops, ordered=False); ops = []
        if ops:
            self.db.bulk_write(ops, ordered=False)

    @staticmethod
    def _running(s, window):
        'centered running mean over `window` samples'
        return s.rolling(max(int(window), 1), center=True, min_periods=1).mean()

    # -------------------------------------------------------------- downloaders
    def _fetch_kpap(self, d0, dn, source='potsdam'):
        rows = []
        if source in ('potsdam', 'potsdam-legacy'):
            txt = _get("https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_Ap_SN_F107_since_1932.txt")
            for ln in txt.splitlines():
                if ln.startswith('#') or len(ln) < 60:
                    continue
                w   = ln.split()                                   # ...Kp1..8  ap1..8  Ap  SN  F10.7obs  F10.7adj  D
                dt  = datetime(int(w[0]), int(w[1]), int(w[2]))
                kp  = [float(x) if float(x) >= 0 else np.nan for x in w[7:15]]
                ap  = [float(x) if float(x) >= 0 else np.nan for x in w[15:23]]
                apd = float(w[23]) if float(w[23]) > 0 else np.nan
                f7  = float(w[25]) if float(w[25]) > 0 else np.nan
                fad = float(w[26]) if float(w[26]) > 0 else np.nan
                sn  = float(w[24]) if float(w[24]) >= 0 else np.nan
                rows.append((dt, kp, ap, apd, np.nanmean(kp), f7, fad, sn))
            if source == 'potsdam-legacy':
                rows = [r[:5] + (np.nan, np.nan, np.nan) for r in rows]     # legacy: no solar flux
        elif source == 'noaa':
            for yr in range(d0.year, dn.year + 1):
                for ln in _get("https://www.ngdc.noaa.gov/stp/GEOMAGNETIC_DATA/INDICES/KP_AP/%i" % yr).splitlines():
                    if len(ln) < 65:
                        continue
                    dt = datetime.strptime(ln[:6].replace(" ", "0"), "%y%m%d")
                    dt = dt.replace(year=dt.year - 100) if dt.year > datetime.now().year else dt
                    kp = [int(ln[12 + 2 * i:14 + 2 * i]) / 10. for i in range(8)]
                    ap = [int(ln[31 + 3 * i:34 + 3 * i]) for i in range(8)]
                    try:    f7 = float(ln[65:70])
                    except Exception: f7 = np.nan
                    rows.append((dt, kp, ap, int(ln[55:58]), np.nanmean(kp), f7, np.nan, np.nan))
        elif source in ('canada', 'lisird'):
            if source == 'canada':
                tab = _get("https://spaceweather.geomag.nrcan.gc.ca/data-donnee/sw-me/sx-5-flux-fluxpeak/daily/fluxtable.txt")
                daily = {}
                for ln in tab.splitlines()[2:]:
                    p = ln.split()
                    if len(p) < 6:
                        continue
                    dt = datetime.strptime(p[0], "%Y%m%d")
                    daily.setdefault(dt, []).append(float(p[5]))
                f107 = {k: np.nanmean([x for x in v if 50 <= x <= 300]) for k, v in daily.items()}
            else:
                f107 = {}
                for yr in range(d0.year, dn.year + 1):
                    for ln in _get("https://lasp.colorado.edu/lisird/latis/dap/noaa_radio_flux.csv?time,f107&time%%3E=%i-01-01&time%%3C%i-01-01" % (yr, yr + 1)).splitlines()[1:]:
                        try:    t, v = ln.split(','); f107[to_datetime(t).to_pydatetime().replace(hour=0)] = float(v)
                        except Exception: pass
            rows = [(dt, [np.nan] * 8, [np.nan] * 8, np.nan, np.nan, v, np.nan, np.nan)
                    for dt, v in f107.items() if _FLR['daily'](d0) <= dt <= dn]
        else:
            raise ValueError("kpap source must be potsdam|potsdam-legacy|noaa|canada|lisird")
        if not rows:
            return
        rows.sort()
        idx    = [dt for dt, *_ in rows]
        f107a  = self._running(Series([r[5] for r in rows], index=idx), 81).values   # 81-day centered mean
        lo, hi = _FLR['daily'](d0) - timedelta(days=100), dn + timedelta(days=100)
        keep   = [i for i, dt in enumerate(idx) if lo <= dt <= hi]
        R      = [rows[i] for i in keep]
        self._write('daily', DataFrame({'datetime': [idx[i] for i in keep], 'f107': [r[5] for r in R],
                                        'f107adj': [r[6] for r in R], 'f107a': [f107a[i] for i in keep],
                                        'ap_daily': [r[3] for r in R], 'kp_daily': [r[4] for r in R],
                                        'sn': [r[7] for r in R]}))
        self._write('h3', concat([DataFrame({'datetime': [r[0] + timedelta(hours=3 * h) for r in R],
                                             'kp': [r[1][h] for r in R], 'ap': [r[2][h] for r in R]})
                                  for h in range(8)]))

    def _fetch_ae(self, d0, dn, source='auto'):
        df = None
        if source in ('auto', 'kyoto'):
            try:    df = self._ae_kyoto(d0, dn)
            except Exception as e: warnings.warn("Kyoto aeasy-cgi failed: %s" % e)
        if source == 'omni':
            warnings.warn("AE from NOAA OMNI 1-min (AE only; AL/AU/AO not populated)")
            df = self._ae_omni(d0, dn)
        elif source in ('auto', 'images'):
            has  = set() if df is None or 'ae' not in df else {t.date() for t in df.index[df['ae'].notna()]}
            days = [d0.date() + timedelta(days=i) for i in range((dn.date() - d0.date()).days + 1)]
            miss = [d for d in days if d not in has] if source == 'auto' else days
            if miss:
                warnings.warn("AE for %i day(s) digitized from Kyoto realtime plots (provisional)" % len(miss))
                for d in miss:
                    try:
                        im = self._ae_images(datetime(d.year, d.month, d.day), datetime(d.year, d.month, d.day))
                        df = im if df is None else df.combine_first(im)
                    except Exception as e:
                        warnings.warn("no realtime AE image for %s (%s)" % (d, e))
        if df is None or df.empty:
            return
        for c in ('ae', 'ao', 'au', 'al'):
            df[c] = df[c] if c in df else np.nan
        df = df[['ae', 'ao', 'au', 'al']].sort_index()
        self._write('min', df.reset_index())
        for w, cad in ((15, 'q15'), (60, 'hourly')):
            g = df.apply(lambda s: self._running(s, w))
            g = g[(g.index.minute % w == 0) if w < 60 else (g.index.minute == 0)]
            self._write(cad, g.add_suffix('_15min' if w == 15 else '_hourly').reset_index())

    def _ae_kyoto(self, d0, dn):
        return self._wdc(_aeasy(d0, dn, 'AE'), slice(21, 23), {'AE': 'ae', 'AO': 'ao', 'AU': 'au', 'AL': 'al'})

    def _fetch_dst(self, d0, dn, *_):
        w = _wdc_lines(_kyoto_month("dst", d0))
        rows = {}
        for ln in w:
            try:
                dt = datetime.strptime(ln[3:5] + ln[5:7] + ln[8:10], "%y%m%d")
                dt = dt.replace(year=dt.year - 100) if dt.year > datetime.now().year else dt
            except Exception:
                continue
            base = tofloat(ln[16:20]) * 100.
            for h in range(24):
                rows[dt + timedelta(hours=h)] = base + tofloat(ln[20 + 4 * h:24 + 4 * h])
        if rows:
            self._write('hourly', DataFrame({'datetime': list(rows), 'dst': list(rows.values())}))

    def _fetch_asy(self, d0, dn, *_):
        df = self._wdc(_aeasy(d0, dn, 'ASY'), slice(21, 24),
                       {'ASYD': 'asyd', 'ASYH': 'asyh', 'SYMD': 'symd', 'SYMH': 'symh'}, suffix_at=18)
        if df is not None and not df.empty:
            self._write('min', df.reset_index())

    # ---- shared parsing -------------------------------------------------------
    def _wdc(self, text, comp_at, want, suffix_at=None):
        'Kyoto WDC hourly-record format -> DataFrame indexed by minute'
        rows = {}
        for ln in text.splitlines():
            if len(ln) < 394:
                continue
            try:
                dt = datetime.strptime(ln[12:18] + " " + ln[19:21], "%y%m%d %H")
                dt = dt.replace(year=dt.year - 100) if dt.year > datetime.now().year else dt
                key = ln[comp_at].strip().upper() + (ln[suffix_at].upper() if suffix_at else "")
            except Exception:
                continue
            col = want.get(key)
            if col is None:
                continue
            for i in range(60):
                rows.setdefault(dt + timedelta(minutes=i), {})[col] = tofloat(ln[34 + 6 * i:40 + 6 * i])
        return DataFrame(rows).T.rename_axis('datetime') if rows else None

    def _ae_omni(self, d0, dn):
        import requests
        p = {'activity': 'retrieve', 'res': 'min', 'spacecraft': 'omni_min', 'vars': '37',
             'start_date': d0.strftime('%Y%m%d'), 'end_date': dn.strftime('%Y%m%d')}
        r = requests.get('https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi', params=p, timeout=120).text
        out = {}
        for ln in r.splitlines():
            w = ln.split()
            if len(w) == 5 and w[0].isdigit():
                dt = datetime(int(w[0]), 1, 1) + timedelta(days=int(w[1]) - 1, hours=int(w[2]), minutes=int(w[3]))
                out[dt] = np.nan if w[4].startswith('9999') else float(w[4])
        return DataFrame({'ae': out}).rename_axis('datetime')

    def _ae_images(self, d0, dn):
        from PIL import Image
        import requests, tempfile
        rows = {}
        dt = datetime(d0.year, d0.month, d0.day)
        while dt <= dn:
            src = dt.strftime("https://wdc.kugi.kyoto-u.ac.jp/ae_realtime/%Y%m/rtae_%Y%m%d.png")
            try:
                raw = requests.get(src, timeout=60)
                if not raw.ok:
                    raise IOError(raw.status_code)
                fp = os.path.join(tempfile.gettempdir(), dt.strftime("rtae_%Y%m%d.png"))
                open(fp, 'wb').write(raw.content)
                for r in _digitize_ae(fp, dt):
                    rows[r[0]] = {'ae': r[1], 'ao': r[2], 'au': r[3], 'al': r[4]}
            except Exception as e:
                warnings.warn("no realtime AE image for %s (%s)" % (dt.date(), e))
            dt += timedelta(days=1)
        return DataFrame(rows).T.rename_axis('datetime') if rows else None


# ---- module-level helpers ---------------------------------------------------
def configure(**kw):
    '''persist connection settings to ~/.config/mpyricalspace/config.json (chmod 600) so
       later DataManager() / manager() calls pick them up without arguments or env vars.
       keys: mongo_uri, mongo_user, mongo_password, mongo_writes, data_dir

       >>> from mpyricalspace.DataManager import configure
       >>> configure(mongo_uri="mongodb://host:27017/", mongo_user="u", mongo_password="p", mongo_writes=True)
    '''
    cfg = _cfg()
    cfg.update({k: v for k, v in kw.items() if v is not None})
    os.makedirs(os.path.dirname(CONFIG), exist_ok=True)
    json.dump(cfg, open(CONFIG, 'w'), indent=2)
    os.chmod(CONFIG, 0o600)
    return CONFIG


DataManager.configure = staticmethod(configure)


def _get(url, timeout=120):
    from urllib.request import urlopen, Request
    with urlopen(Request(url, headers={'User-Agent': 'mpyricalspace'}), timeout=timeout) as r:
        return r.read().decode('utf-8', 'replace')


_IRI_INDEX_URL = 'https://chain-new.chain-project.net/echaim_downloads/'


def fetch_iri_indices(*dests):
    'download the current ig_rz.dat + apf107.dat (CHAIN / ECHAIM mirror) into each dest directory'
    for name in ('ig_rz.dat', 'apf107.dat'):
        txt = _get(_IRI_INDEX_URL + name, timeout=90)
        for d in dests:
            tmp = os.path.join(d, name + '.tmp')
            open(tmp, 'w').write(txt)
            os.replace(tmp, os.path.join(d, name))            # atomic; breaks any hard-link


def _kyoto_month(kind, d0):
    ym = d0.strftime("%Y%m")
    fn = "%s%s%02i.for.request" % (kind, d0.strftime("%y"), d0.month)
    for stage in ("final", "provisional", "realtime"):
        try:    return _get("https://wdc.kugi.kyoto-u.ac.jp/%s_%s/%s/%s" % (kind, stage, ym, fn))
        except Exception: continue
    raise IOError("no Kyoto %s file for %s" % (kind, ym))


_wdc_lines = lambda t: [ln for ln in t.splitlines() if ln.strip()]


def _aeasy(d0, dn, output):
    'Kyoto aeasy-cgi -> WDC text for a whole month (output: AE | ASY)'
    import requests
    from collections import OrderedDict
    dur = (datetime(d0.year + (d0.month == 12), d0.month % 12 + 1, 1) - datetime(d0.year, d0.month, 1)).days - 1
    q = OrderedDict([('Tens', '%03i' % (d0.year // 10)), ('Year', '%1i' % (d0.year % 10)), ('Month', '%02i' % d0.month),
                     ('Day_Tens', '0'), ('Days', '1'), ('Hour', '00'), ('min', '00'),
                     ('Dur_Day_Tens', '%02i' % (dur // 10)), ('Dur_Day', '%i' % (dur % 10)),
                     ('Dur_Hour', '23'), ('Dur_Min', '59'), ('Image Type', 'GIF'), ('COLOR', 'COLOR'),
                     ('AE Sensitivity', '0'), ('ASY/SYM  Sensitivity', '0'),
                     ('Output', output), ('Out format', 'WDC'), ('Email', 'mpyricalspace@example.com')])
    t = requests.get("https://wdc.kugi.kyoto-u.ac.jp/cgi-bin/aeasy-cgi", params=q, timeout=180).text
    if '<html' in t.lower() or '<HTML' in t:
        raise IOError("aeasy-cgi returned HTML (no data for %s)" % d0.strftime('%Y-%m'))
    return t


def _digitize_ae(path, dt):
    'read AE/AO/AU/AL 1-min values off a Kyoto rtae_YYYYMMDD.png'
    from PIL import Image
    img = np.asarray(Image.open(path).convert('L'))          # 0=black trace, 255=white background
    lower = _clean_pixels(img[254:396, 81:649])
    upper = _clean_pixels(img[27:198, 81:649])

    def scale(top, bottom, delta, h, w):
        rdx = np.round(np.arange(0, h, (h - 1) / 5.)).astype(int)
        aes, d = [], 0.
        for i in range(rdx[-1]):
            if i in rdx: d = delta / (rdx[list(rdx).index(i) + 1] - rdx[list(rdx).index(i)])
            aes.append(top if i == 0 else aes[-1] - d)
        aes.append(bottom)
        cdx = np.round(np.arange(0, w, (w - 1) / 24.)).astype(int)
        for k, v in ((4, 95), (7, 166), (-10, 355), (-5, 473), (-2, 544)):
            cdx[k] = v
        dts, dd = [], 0.
        for i in range(cdx[-1]):
            if i in cdx: dd = 60. / (cdx[list(cdx).index(i) + 1] - cdx[list(cdx).index(i)])
            dts.append(dt if i == 0 else dts[-1] + timedelta(minutes=dd))
        dts.append(dt + timedelta(hours=24))
        return dts, {i: a for i, a in enumerate(aes)}

    dts, lut1 = scale(2000., -500., 500., 142, 568)
    _,   lut2 = scale(1000., -2000., 500., 171, 568)
    vals = []
    for j in range(lower.shape[1]):
        c1 = np.where(np.isfinite(lower[:, j]))[0]
        c2 = np.where(np.isfinite(upper[:, j]))[0]
        ae, ao = (lut1[c1[0]], lut1[c1[-1]]) if len(c1) else (np.nan, np.nan)
        au, al = (lut2[c2[0]], lut2[c2[-1]]) if len(c2) else (np.nan, np.nan)
        hrs = dts[j].hour + dts[j].minute / 60. + dts[j].second / 3600.
        vals.append([hrs + (24. if hrs <= 0 and dts[j] > dts[0] else 0.), ae, ao, au, al])
    v = np.array(vals)
    mins = [h + m / 60. for h in range(24) for m in range(60)]
    out = []
    cols = [np.interp(mins, v[np.isfinite(v[:, k]), 0], v[np.isfinite(v[:, k]), k]) for k in range(1, 5)]
    for i, (h, mm) in enumerate([(h, m) for h in range(24) for m in range(60)]):
        out.append([dt + timedelta(hours=h, minutes=mm)] + [c[i] for c in cols])
    return out


def _clean_pixels(data):
    'keep the dark plot trace, drop background + thin noise'
    d = np.where(data < 128, data.astype(float), np.nan)           # trace = dark pixels
    for _ in range(3):
        nan = ~np.isfinite(d)
        drop = np.zeros_like(nan)
        drop[:, 1:-1] |= nan[:, :-2] & nan[:, 2:]                  # nan on both sides horizontally
        drop[1:-1, :] |= nan[:-2, :] & nan[2:, :]                  # ... or vertically
        d[np.isfinite(d) & drop] = np.nan
    return d


_DEFAULT = None
def manager():
    'process-wide default DataManager (used by Predictor.py)'
    global _DEFAULT
    if _DEFAULT is None:
        _DEFAULT = DataManager()
    return _DEFAULT
