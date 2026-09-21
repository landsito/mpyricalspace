'''
Empirical: a thin facade over mpyricalspace.models.

    obj = Empirical()
    obj.set_period(d0, dn, freq=timedelta(minutes=15))   # or obj.set_time([...])
    obj.set_location(lat, lon, alt)                      # needed by hwm / msis / igrf
    ds = obj.run_rocsat()                                # every run_* returns a fresh xr.Dataset

run_* are pure: they read obj.time / obj.lat|lon|alt and return a Dataset, they keep
no state on obj. Combine results with xarray.merge. The models themselves live in
models.py and can be called directly.
'''
import numpy as np
import xarray as xr
from datetime import datetime, timedelta

from mpyricalspace import models
from mpyricalspace.DataManager import manager


class Empirical(object):

    def set_period(self, d0=datetime(2010, 6, 20), dn=datetime(2010, 6, 21), freq=timedelta(minutes=15.)):
        self.d0, self.dn, self.freq = d0, dn, freq
        self.time = np.arange(d0, dn + freq, freq)

    def set_time(self, dts):
        self.time = np.atleast_1d(np.asarray(dts, dtype='datetime64[s]'))
        t = self.time.astype('datetime64[s]').astype(datetime)
        self.d0, self.dn = t.min(), t.max()
        self.freq = timedelta(0) if self.time.size < 2 else (t[1] - t[0])

    def set_location(self, lat, lon, alt):
        self.lat, self.lon, self.alt = lat, lon, alt

    # ---------------------------------------------------------------- facade
    def run_scherliessfejer(self, hour_resolution=False, F107=None, indices=None, SLT=None,
                            half_resolution=False):
        '''Scherliess-Fejer quiet + storm vertical drift, evaluated at every obj.time.
           `half_resolution` is deprecated and ignored.'''
        return models.scherliess_fejer(self.time, lon=getattr(self, 'lon', None),
                                       hour_resolution=hour_resolution, F107=F107, indices=indices,
                                       SLT=SLT, half_resolution=half_resolution)

    def run_jvdm1(self, f107=None, f107a=None, slt=None, doy=None):
        return models.jvdm1_drift(self.time, f107=f107, f107a=f107a, slt=slt, doy=doy)

    def run_rocsat(self, f107s=None, lons=None, doys=None, slts=None):
        return models.rocsat_drift(self.time, lons if lons is not None else getattr(self, 'lon', None),
                                   f107s=f107s, doys=doys, slts=slts)

    def run_eej(self, version=2, lon=None, flux=None, slts=None, doys=None, lunars=None,
                model="champ"):
        '''Alken Equatorial Electrojet climatology -- version = 2 (default) | 1.
           model (v2 only): "champ" (default) | "oersted" | "sac-c" -- the satellite
           the coefficients were fit from (Alken & Maus 2007).'''
        return models.eej(self.time, lon if lon is not None else getattr(self, 'lon', None),
                          flux=flux, slts=slts, doys=doys, lunars=lunars, version=version,
                          model=model)

    def run_eef(self, lon=None, flux=None, slts=None, doys=None, lunars=None):
        '''Alken Equatorial Electric Field climatology.'''
        return models.eef(self.time, lon if lon is not None else getattr(self, 'lon', None),
                          flux=flux, slts=slts, doys=doys, lunars=lunars)

    def run_hwm(self, version=2014, ap=None, lon=None, lat=None, alt=None, ut=None, doy=None):
        '''Horizontal Wind Model -- version = 2014 (default) | 2007 | 1993.'''
        return models.hwm(self.time, self.lat if lat is None else lat, self.lon if lon is None else lon,
                          self.alt if alt is None else alt, ap=ap, version=version, ut=ut, doy=doy)

    def run_hltwim(self, kp=None, lon=None, lat=None, ut=None, doy=None):
        '''High-Latitude Thermospheric Wind Model (Dhadly et al. 2019); |MLAT|>40 only.'''
        return models.hltwim(self.time, self.lat if lat is None else lat,
                             self.lon if lon is None else lon, kp=kp, ut=ut, doy=doy)

    def run_msis(self, ap=None, lon=None, lat=None, alt=None, f107s=None,
                nativelypackaged_indices=False, rho_m3=True):
        '''nativelypackaged_indices=False (default): F10.7/F10.7a/ap from DataManager.
           nativelypackaged_indices=True: from pymsis's own bundled/cached index file
           (pymsis.utils.get_f107_ap, CelesTrak) instead. See models.msis.'''
        return models.msis(self.time, self.lat if lat is None else lat, self.lon if lon is None else lon,
                           self.alt if alt is None else alt, f107s=f107s, ap=ap,
                           nativelypackaged_indices=nativelypackaged_indices, rho_m3=rho_m3)

    def run_weimer05(self, mlat, mlt, by=None, bz=None, vsw=None, nsw=None, tilt=None, res='5min', avg=20, lag=0):
        '''
        Weimer (2005) high-latitude electric potential [kV] and field-aligned current [uA/m^2, + = downward]
        at obj.time, on AACGM latitudes `mlat` [deg; negative = southern hemisphere] and magnetic local times
        `mlt` [h] -- see models.weimer05 for the coordinates, the (time, mlat, mlt) cube / aligned-samples rule,
        and NaN handling.

        AVERAGING DEFAULT: drivers not passed in (by, bz [nT, GSM], vsw [km/s], nsw [cm^-3]) come from the
        DataManager as the MEAN OF THE PREVIOUS 20 MIN of the 5-min OMNI series (res='5min', avg=20), the way
        Weimer (2005b) drives the model -- not as instantaneous values. Pass avg=None, res='1min' for the
        instantaneous 1-min value, or avg=45, lag=10 for the recipe the coefficients were fitted with
        (Weimer 2005a). tilt=None computes the dipole tilt from obj.time.
        '''
        return models.weimer05(self.time, mlat, mlt, by=by, bz=bz, vsw=vsw, nsw=nsw, tilt=tilt,
                               res=res, avg=avg, lag=lag)

    def run_igrf(self, version='latest', lats=None, lons=None, alts=None):
        return models.igrf(self.time, getattr(self, 'lat', None) if lats is None else lats,
                           getattr(self, 'lon', None) if lons is None else lons,
                           getattr(self, 'alt', None) if alts is None else alts, version=version)

    def run_iri(self, lats=None, lons=None, alts=None, quiet=True, F107=None, F107a=None,
                nativelypackaged_indices=False, NmF2=None, hmF2=None,
                version=2026, compute_Ne=True, compute_Te_Ti=True, compute_Ni=True, rho_m3=True):
        '''nativelypackaged_indices=False (default): F10.7/F10.7a from DataManager.
           nativelypackaged_indices=True: F10.7/F10.7a from IRI's own bundled apf107.dat.
           Rz12/IG12 always come from IRI's own ig_rz.dat either way. See models.iri.'''
        return models.iri(self.time, lats, lons, alts, quiet=quiet, F107=F107, F107a=F107a,
                          nativelypackaged_indices=nativelypackaged_indices, NmF2=NmF2, hmF2=hmF2,
                          version=version, compute_Ne=compute_Ne, compute_Te_Ti=compute_Te_Ti,
                          compute_Ni=compute_Ni, rho_m3=rho_m3)

    def run_ppeefm1(self, data=None, nativelypackaged_code=False, ief=None,
                    nativelypackaged_indices=False):
        '''Manoj & Maus prompt-penetration equatorial electric field, as a vertical drift.

           nativelypackaged_code=False (default) hits the Manoj & Maus PPEF web service;
           nativelypackaged_code=True runs the same model (RTEEF) locally, from the
           bundled ACE data (2001-2007 only). Pass a (datetime, qef, tef, ppef) array,
           or a dict {'lla': (lat, lon, alt), 'data': array}, to use your own.

           ief=<1-D array of synthetic IEF Ey [mV/m]> feeds a made-up solar-wind
           electric field straight into the model's TF.COF prompt-penetration filter
           (it overrides the internal ACE buffer). 5-min cadence, length
           (dn - d0)/300 s + 13 -- one hour of priming before obj.time[0], through
           obj.time[-1]. Returns the raw model output in mV/m (qef, ppef, tef), not
           the drift-scaled fields; qef is 0 in this mode.

           nativelypackaged_indices (nativelypackaged_code=True or ief= only): the F10.7
           the climatology (qef) is built from -- False (default) DataManager, True the
           bundled SPIDR F107.txt; either way it only sets the baseline, not the
           prompt part. No effect when nativelypackaged_code=False.'''
        if isinstance(data, dict):
            self.set_location(*data['lla'])
            data = data['data']
        return models.manoj_maus(self.time, self.d0, self.dn, self.freq,
                                 lon=getattr(self, 'lon', -77.), data=data,
                                 nativelypackaged_code=nativelypackaged_code, ief=ief,
                                 nativelypackaged_indices=nativelypackaged_indices)

    # -------------------------------------------------- multi-model (survey.py)
    def run_track(self, lats, lons, alts=None, models=None, equator_deg=20.0, **opts):
        '''Every applicable model along a trajectory sampled at obj.time -- one merged
           xr.Dataset. lats/lons/alts are 1-D, the same length as obj.time (a scalar
           broadcasts). models=None (default) runs them all; pass a list of names
           (survey.model_domains() for the allowed names). See
           mpyricalspace.survey.run_track.'''
        from mpyricalspace import survey
        return survey.run_track(self.time, lats, lons, alts, models=models,
                                equator_deg=equator_deg, **opts)

    def run_grid(self, lats=None, lons=None, alts=None, models=None, equator_deg=20.0, **opts):
        '''Every applicable model on the meshgrid obj.time x lats x lons x alts.
           See mpyricalspace.survey.run_grid.'''
        from mpyricalspace import survey
        return survey.run_grid(times=self.time, lats=lats, lons=lons, alts=alts,
                               models=models, equator_deg=equator_deg, **opts)


if __name__ == "__main__":
    obj = Empirical()
    obj.set_period(datetime(2024, 5, 11, 20), datetime(2024, 5, 12, 12), timedelta(minutes=5))
    obj.set_location(*models.JRO)
    print(obj.run_scherliessfejer())
