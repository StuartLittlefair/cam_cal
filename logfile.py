import numpy as np
import numpy.ma as ma
from hipercam import hlog
from astropy.io import fits
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, name_resolve
import os
import re
from tkinter.filedialog import askopenfilename, askopenfilenames


class Logfile:

    def __init__(self, file, instrument, tel_location, target_coords=None, verbose=True):
        self.logfile = file
        self.instrument = instrument
        self.tel_location = tel_location

        self.run = self.getRun()
        self.path, self.fname = self.getPath()
        self.target, self.filters, self.target_coords = self.getTarget(target_coords, verbose)
        self.logf = hlog.Hlog.rascii(self.logfile)
        self.apnames = self.logf.apnames.copy()
        if verbose:
            print('Target = {}'.format(self.target))
            print('Run = {}'.format(self.run))
            print('Filters = {}'.format(self.filters))


    def getPath(self):
        path = os.path.normpath(self.logfile)
        directory_tree = path.split(os.sep)
        mainpath = os.path.split(self.logfile)[0]
        fname = directory_tree[-1]
        return mainpath, fname


    def getRun(self):
        log_params = dict()
        fptr = open(self.logfile)
        for line in fptr:
            eq = line.find('=')
            if eq > -1:
                log_params[line[:eq].strip('#').strip()] = line[eq+1:].split('#')[0].strip()
        fptr.close()
        return log_params['run']


    def getTarget(self, target_coords, verbose=True):
        f = fits.open(os.path.join(self.path, self.run + '.hcm'))
        if self.instrument == 'ultracam':
            target = f[0].header['TARGET']
            filters = re.findall(r"\w?\s?([ugriz])'", f[0].header['FILTERS'])
            if target_coords:
                coords = target_coords
            else:
                coords = self.getCoords(target, verbose=verbose)
        elif self.instrument == 'hipercam':
            target = f[0].header['OBJECT']
            filters = re.findall(r"([ugriz])s,?", f[0].header['FILTERS'])
            coords_str = ' '.join([f[0].header['RA'], f[0].header['DEC']])
            coords = SkyCoord(coords_str, unit=(u.hourangle, u.deg))
        f.close()
        return target, filters, coords
    

    def manualInput(self):
        obj_simbad = input('Enter Simbad object name '
                        '(to enter coordinates press <ENTER>): ')
        if obj_simbad != '':
            target_coords = SkyCoord.from_name(obj_simbad, parse=True)
            print('Found object\n')
        else:
            coords = input('Enter coordinates: ')
            target_coords = SkyCoord(coords, unit=(u.hourangle, u.deg))
        return target_coords


    def getCoords(self, target=None, verbose=True):
             if target:
                 try:
                     target_coords = SkyCoord.from_name(target, parse=True)
                 except name_resolve.NameResolveError:
                     print(f"Can't find {target} in Simbad, enter object name manually")
                     target_coords = self.manualInput()
             else:
                 target_coords = self.manualInput()
             return target_coords


    def openData(self, ccd, ap, save=False, mask=True):
        target = self.target.replace(' ', '_')
        filters = self.filters[::-1]
        ccd = str(ccd)
        ap = str(ap)
        if ccd not in self.apnames.keys():
            raise ValueError('{} not a valid CCD'.format(ccd))
        if ap not in self.apnames[ccd]:
            raise ValueError('{} not a valid aperture'.format(ap))

        data = self.logf.tseries(ccd, ap)
        obstime = Time(data.t,format='mjd', scale='utc',location=self.tel_location)
        data.t = obstime.tdb.value
        exp = self.logf[ccd]['Exptim']/86400
        weights = np.ones(len(data.t))
        m = data.get_mask()
        zero_flux_mask = ma.getmask(ma.masked_not_equal(data.y, 0))
        data_mask = np.logical_and(~m, zero_flux_mask)
        if mask:
            out = np.column_stack([data.t[data_mask],
                                   exp[data_mask],
                                   data.y[data_mask],
                                   data.ye[data_mask],
                                   weights[data_mask],
                                   weights[data_mask]])
        else:
            out = np.column_stack([data.t, exp, data.y, data.ye, weights, weights])
        return out, data_mask
    

    def barycorr(self, t):
        obstimes = Time(t, format='mjd', scale='tdb', location=self.tel_location)
        bary_corr = obstimes.light_travel_time(self.target_coords)
        return bary_corr