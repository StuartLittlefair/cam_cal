import numpy as np
from hipercam import hlog
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy.ma as ma
from scipy.optimize import curve_fit, leastsq
from scipy.interpolate import interp1d
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, name_resolve
import astropy.units as u
from astropy.time import Time
from tkinter.filedialog import askopenfilename, askopenfilenames
from tkinter import Tk
import os
import re
import warnings
from mergedeep import merge
import heapq


warnings.filterwarnings('error')


def linear(x, m, c):
    return m*x + c


def flux_to_mag(flux, error):
    mag = -2.5*np.log10(flux)
    mag_err = (2.5 / np.log(10)) * (error/flux)
    return mag, mag_err


def mag_to_flux(mag):
    flux = 10**(-0.4 * mag)
    return flux


def magAB_to_flux(mag, mag_err=0):
    ABflux = 10**(-0.4 * (mag-8.90)) * u.Jy
    ABflux_err = mag_err * (np.log(10) / 2.5) * ABflux
    return ABflux, ABflux_err

def flux_to_ABmag(flux, flux_err=0):
    ABmag = -2.5 * np.log10(flux.to_value(u.Jy)) + 8.90
    ABmagErr = (flux_err.value/flux.value) * (2.5 / np.log(10))
    return ABmag, ABmagErr

def clip(fit, x, y, sigma):
    std = np.std(y - fit(x))
    med = np.median(y - fit(x))
    mask = ma.getmask(ma.masked_inside(y - fit(x), med-sigma*std, med+sigma*std))
    return mask


def polysigclip(x, y, porder, sigma, iters):
    mask = np.ones(x.shape, dtype=bool)
    for _ in range(iters):
        pfit = np.polyfit(x[mask], y[mask], porder)
        p = np.poly1d(pfit)
        std = np.std(y[mask] - p(x[mask]))
        mask = ma.getmask(ma.masked_inside(y - p(x), -sigma * std, sigma * std))
    return mask, p


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
        if save and get_target:
            if not os.path.isdir(os.path.join(self.path, 'reduced')):
                os.makedirs(os.path.join(self.path, 'reduced'))
            fname = '{}_{}_{}_ap{}.dat'.format(target, self.run, filters[int(ccd)-1], ap)
            fname = os.path.join(self.path, 'reduced', fname)
            np.savetxt(fname, out, fmt='%11.9f %9.4e %9.4e %9.4e %1.0f %1.0f')
        return out, data_mask
    

    def barycorr(self, t):
        obstimes = Time(t, format='mjd', scale='tdb', location=self.tel_location)
        bary_corr = obstimes.light_travel_time(self.target_coords)
        return bary_corr
        


class Observation:
    """Class for reducing and flux calibrating ULTRACAM and HiPERCAM
       eclipse light curves."""


    def __init__(self, instrument='ultracam', tel_location=None):
        root = Tk()
        root.withdraw()
        self.instrument = instrument
        # set instrument specific variables
        if self.instrument=='ultracam':
            self.tel_location = EarthLocation.of_site('La Silla Observatory')
            self.filt2ccd = dict(u='3', g='2', r='1', i='1', z='1')
            self.rootDataDir = '~/Astro/Data/photometry/ultracam'
        elif self.instrument=='hipercam':
            self.tel_location = EarthLocation.of_site('Roque de los Muchachos')
            self.filt2ccd = dict(u='1', g='2', r='3', i='4', z='5')
            self.rootDataDir = '~/Astro/Data/photometry/hipercam'
        else: raise ValueError('{} is not a valid instrument'.format(self.instrument))
        
        if tel_location:
            if tel_location in EarthLocation.get_site_names():
                self.tel_location = EarthLocation.of_site(tel_location)
            elif isinstance(tel_location, EarthLocation):
                self.tel_location = tel_location
            else:
                raise ValueError('"{}" not a valid EarthLocation site name or '
                                 'an EarthLocation object'.format(tel_location))
        # set extinction coeffs and instrumental zeropoints
        self.set_default()
        #initialise observation dict
        self.observations = dict.fromkeys(['science', 'std', 'atm'])

    
    def set_default(self):
        """Sets default atmospheric extinctions and instrument zeropoints"""

        self.atm_extinction = dict(mean=dict(u=0.48, g=0.24, r=0.18, i=0.15, z=0.10),
                                   err=dict(u=0.05, g=0.05, r=0.05, i=0.05, z=0.05))
        if self.instrument=='ultracam':
            self.zeropoint = dict(mean=dict(u=25.09, g=26.70, r=26.30, i=25.90, z=25.30),
                                  err=dict(u=0.05, g=0.05, r=0.05, i=0.05, z=0.05),
                                  airmass=dict(airmass=0))
        elif self.instrument=='hipercam':
            self.zeropoint = dict(mean=dict(u=28.15, g=29.22, r=28.78, i=28.43, z=27.94),
                                  err=dict(u=0.05, g=0.05, r=0.05, i=0.05, z=0.05),
                                  airmass=dict(airmass=0))


    def add_observation(self, name=None, logfiles=None, obs_type='science', cal_mags=None):
        """Adds logfiles and any other key info about a target/observation to
           the observations dictionary. Observation type is specified as 'science'
           for science target data, 'std' for flux standard observations, and 'atm'
           for long runs for measuring atmospheric extinction."""

        if not obs_type:
            obs_type = input("'science', 'std', or 'atm'? ")
        if not name:
            name = input('Enter name for observation: ')
        if not logfiles:
            logfiles = list(askopenfilenames(title='Select observation logfiles',
                                             initialdir=self.rootDataDir,
                                             filetypes=[("hcam logfile", ".log")]))
        if obs_type=='std' and not cal_mags:
            raise ValueError('If observations are of a flux standard then '
                             'calibrated magnitudes must be supplied')
        obs1, obs2, obs3 = dict(), dict(), dict()
        obs1['name'], obs1['logfiles'] = name, logfiles
        if obs_type == 'std':
            obs1['cal_mags'] = cal_mags
        obs2[name] = obs1
        obs3[obs_type] = obs2
        self.observations = merge(self.observations, obs3)
    

    def remove_observation(self, name, obs_type=None):
        """removes observation entry from observation dictionary."""

        if not obs_type:
            for obstype in list(self.observations.keys()):
                if self.observations[obstype]:
                    del self.observations[obstype][name]
        else:
            del self.observations[obs_type][name]
    

    def clear(self):
        self.__init__(instrument=self.instrument)


    def airmass(self, times, target_coords):
        """Returns target airmass given an array of astropy times and the target position"""

        frame = AltAz(obstime=times, location=self.tel_location)
        altazs = target_coords.transform_to(frame)
        airmasses = altazs.secz
        return airmasses


    def get_atm_ex(self, plot=True):
        """calulates atmospheric extinction from all apertures included in 'atm' logfiles.
           If two logfiles are of the same field then it will stich the runs together and fit as one.
           Apertures must be totally consistent between runs."""

        if 'atm' not in self.observations.keys():
            raise ValueError('No atmospheric extinction data added.')
        atm_targets = dict()
        temp_extinctions = dict(mean=dict(), err=dict())
        # get logfiles of all added 'atm' observations
        logfiles = [self.observations['atm'][name]['logfiles']
                    for name in self.observations['atm'].keys()]
        # flatten list of lists
        logfiles = [item for sublist in logfiles for item in sublist]

        for file in logfiles:
            log = Logfile(file, self.instrument, self.tel_location, verbose=False)
            filters = log.filters
            if log.target not in atm_targets.keys():
                atm_targets[log.target] = [file]
            else:
                atm_targets[log.target].append(file)

        for filt in filters:
            ap_info = []
            apertures = []

            for target in atm_targets.keys():
                log = Logfile(atm_targets[target][0], self.instrument,
                              self.tel_location, verbose=False)

                for ap in log.apnames[self.filt2ccd[filt]]:
                    data = np.empty((0,4))

                    for lfile in atm_targets[target]:
                        # stitch multiple observations of same field together
                        log = Logfile(lfile, self.instrument, self.tel_location, log.target_coords, verbose=False)
                        data_new = log.openData(self.filt2ccd[filt], ap,
                                                save=False,
                                                mask=True)[0][:, [0, 2, 3]]
                        t = Time(data_new[:,0], format='mjd', scale='tdb')
                        # add airmass column
                        data_am = np.column_stack((self.airmass(t, log.target_coords), data_new))
                        data = np.vstack((data, data_am))
                        # remove negative flux measurements
                        data = data[data[:,2] > 0]
                    ap_info.append((target, ap))
                    # add stitched data to apertures list
                    apertures.append(data)
            p, errs = self.__fit_atm_ext__(filt, apertures)
            print('{}-band extinction: {:.3f} +- {:.3f}'.format(filt, p[0], errs[0]))
            if plot:
                _, ax = plt.subplots()
                for idx, apdata in enumerate(apertures):
                    ax.scatter(apdata[:,0], flux_to_mag(apdata[:,2], apdata[:,3])[0],
                               label="ap{}, {}".format(*ap_info[idx]))
                    ax.plot(apdata[:,0], linear(apdata[:,0], p[0], p[idx+1]), 'k--')
                plt.show()
            temp_extinctions['mean'][filt] = p[0]
            temp_extinctions['err'][filt] = np.max([errs[0], 0.01])
        write = input("Set fitted atmospheric extinction values [y/n]? ")
        if not write or write=='y':
            # overwrite default extinctions with fitted measurements
            for key, value in temp_extinctions.items():
                self.atm_extinction[key] = value
            print('Atmospheric extinction set.\n')
        else:
            print('Keeping default atmospheric extinction.\n')


    def __atm_chisq__(self, pars, apertures):
        """Computes combined chi (not squared) for every aperture"""

        k, *c = pars
        diff = np.array([])
        for i, data in enumerate(apertures):
            # compute chi residuals for each aperture
            mag, mag_err = flux_to_mag(data[:,2], data[:,3])
            airmass = data[:,0]
            diff_new = np.abs((mag - linear(airmass, k, c[i])) / mag_err)
            diff = np.hstack((diff, diff_new))
        return diff
    

    def __red_chisq__(self, pars, apertures):
        """Reduced chi-squared given data and fit parameters"""

        chisq = np.sum(self.__atm_chisq__(pars, apertures)**2)
        ndata = 0
        for data in apertures:
            ndata += data.shape[0]
        ndf = ndata-len(pars)
        red_chisq = chisq / ndf
        return red_chisq.value


    def __fit_atm_ext__(self, filt, apertures):
        """Fits atmospheric extinction coefficient to data in given apertures."""
        
        # add offset parameter for each aperture
        x0 = [self.atm_extinction['mean'][filt]] + [0 for i in apertures]
        p, cov, _, _, _ = leastsq(self.__atm_chisq__, x0,
                                  full_output=True, args=(apertures))
        red_chisq = self.__red_chisq__(p, apertures)
        # scale variances with reduced chisqr and get std_err for fitted params
        errs = np.sqrt(np.diag(cov * red_chisq))
        return p, errs


    def get_zeropoint(self):
        """Calculate magnitude of star that would give 1 count/sec at the
           instrument if no atmosphere. Must have added a 'std' observation
           with calibrated magnitudes in the filter system used."""

        # Need to check multiple logfiles for same target at different airmasses
        # check zip(target_names, logfiles) will always be the same length

        if 'std' not in self.observations.keys():
            raise ValueError('No atmospheric extinction data added.')
        zp_dict = dict(mean=dict(), err=dict())
        # get logfiles of all added 'std' observations
        logfiles = [self.observations['std'][name]['logfiles']
                    for name in self.observations['std'].keys()]
        # flatten list of lists
        logfiles = [item for sublist in logfiles for item in sublist]
        target_names = [name for name in self.observations['std'].keys()]
        for name, file in zip(target_names, logfiles):
            log = Logfile(file, self.instrument, self.tel_location)
            for filt in log.filters:
                data = log.openData(self.filt2ccd[filt],
                                    ap='1')[0][:, [0, 1, 2, 3]]
                t, te, y, ye = np.split(data, 4, axis=1)
                flux = y/(te*86400)
                flux_err = ye/(te*86400)
                flux_std_err = np.mean(flux_err) / np.sqrt(len(flux_err))
                t = Time(t, format='mjd', scale='tdb')
                airmass = self.airmass(t, log.target_coords)
                mag_i, mag_i_err = flux_to_mag(flux, flux_std_err)
                mag_i0 = mag_i - (self.atm_extinction['mean'][filt] * airmass)
                mag_i0_err = (mag_i_err**2
                              + (self.atm_extinction['err'][filt]
                                 * airmass)**2
                             )**0.5
                zp = self.observations['std'][name]['cal_mags']['mean'][filt] - np.mean(mag_i0)
                zp_err = (np.mean(mag_i0_err)**2
                          + self.observations['std'][name]['cal_mags']['err'][filt]**2)**0.5
                zp_dict['mean'][filt] = zp.value
                zp_dict['err'][filt] = zp_err.value
                zp_dict['airmass'] = np.mean(airmass)
                print('{}-band zeropoint: {:.3f} +- {:.3f} (airmass = {:.2f})'
                      ''.format(filt, zp, zp_err, np.mean(airmass)))
        write = input("Set calculated zeropoints [y/n]? ")
        if not write or write=='y':
            for key, value in zp_dict.items():
                self.zeropoint[key] = value
            print('Zeropoints set.\n')
        else:
            print('Keeping default zeropoints.\n')

    
    def __flux_cal_err__(self, airmass_t, mag_i_err, filt):
        mstd_mstdi_var = (self.zeropoint['err'][filt]**2
                          - (self.zeropoint['airmass']
                             * self.atm_extinction['err'][filt])**2)
        am_diff = self.zeropoint['airmass'] - airmass_t
        airmass_var = (self.atm_extinction['err'][filt] * am_diff)**2
        mag_t_err = np.sqrt( mstd_mstdi_var + airmass_var +  mag_i_err**2 )
        return mag_t_err



    def __cal_comp__(self, data, filt, coords):
        """Calibrate the mean flux of the comparison and it's uncertainty using
           the instrumental zeropoint and atmospheric extinction"""
        # TODO: implement outputting/storing comparison magnitude.
        t, te, y, ye, _, _ = np.split(data, 6, axis=1)
        flux = y/(te*86400)
        flux_err = ye/(te*86400)
        t = Time(t, format='mjd', scale='tdb')
        airmass = self.airmass(t, coords)
        mag_i, mag_i_err = flux_to_mag(flux, flux_err)
        mag_i0 = mag_i - (self.atm_extinction['mean'][filt] * airmass)
        mag_comp = mag_i0 + self.zeropoint['mean'][filt]
        mag_comp_err = self.__flux_cal_err__(airmass, mag_i_err, filt)
        snr_comp = np.median(flux/flux_err)
        comp_flux, comp_flux_err = magAB_to_flux(np.mean(mag_comp), np.mean(mag_comp_err))
        return comp_flux, comp_flux_err, snr_comp, np.mean(airmass)
    

    def __get_eclipse__(self, time, flux, width=1):
        """Chop out width x eclipse_width either side of the middle of eclipse.
           width=1 gives the eclipse plus half the eclipse width either side."""

        # half_flux = (np.max(flux) + np.min(flux)) / 2
        half_flux = (np.median(heapq.nlargest(10, flux)) + np.median(heapq.nsmallest(10, flux))) / 2
        a1 = flux.copy()
        idx1 = np.nanargmin(np.abs(a1 - half_flux))
        a1[idx1-4:idx1+5] = np.nan
        idx2 = np.nanargmin(np.abs(a1 - half_flux))
        print(half_flux, idx1, idx2)
        fig, ax = plt.subplots()
        ax.scatter(time, flux)
        ax.scatter(time[idx1-4:idx1+5], flux[idx1-4:idx1+5])
        ax.scatter(time[idx2-4:idx2+5], flux[idx2-4:idx2+5])
        ax.scatter(time[idx1], flux[idx1])
        ax.scatter(time[idx2], flux[idx2])
        plt.show()
        interp1 = interp1d(flux[idx1-4:idx1+5], time[idx1-4:idx1+5])
        interp2 = interp1d(flux[idx2-4:idx2+5], time[idx2-4:idx2+5])
        t1 = interp1(half_flux)
        t2 = interp2(half_flux)
        eclipsewidth = np.abs(t2-t1)
        mideclipse = (t1 + t2) / 2
        print(mideclipse)
        start = mideclipse - (width * eclipsewidth)
        end = mideclipse + (width * eclipsewidth)
        return start, end
    
    def calibrate_science(self, target_name, eclipse=True):
        """Flux calibrate the selected science target using the calibrated comparison stars."""
    # TODO: implement supplying comparison star magnitude so many runs can be calibrated using one good run.

        log = Logfile(self.observations['science'][target_name]['logfiles'][0],
                      self.instrument, self.tel_location)
        target = log.target.replace(' ', '_')
        ccd = self.filt2ccd['g']
        target_data_orig, target_mask = log.openData(ccd, '1', save=False,
                                                     mask=False)
        comp_data, comp_mask = log.openData(ccd, '2', save=False,
                                            mask=False)
        master_mask = np.logical_and(target_mask, comp_mask)
        target_data = target_data_orig[master_mask]
        comp_data = comp_data[master_mask]
        t_t, t_te, t_y, t_ye, _, _ = ((map(np.squeeze, np.split(target_data, 6, axis=1))))
        _, _, c_y, c_ye, _, _ = ((map(np.squeeze, np.split(comp_data, 6, axis=1))))
        diffFlux = t_y / c_y

        if eclipse:
            start, end = self.__get_eclipse__(t_t, diffFlux, width=1)
        else:
            start, end = t_t[0]-0.0001, t_t[-1]+0.0001 

        for filt in log.filters:
            out_dict = dict(data=dict(), fname=dict())
            ccd = self.filt2ccd[filt]
            print('\n{}-band (CCD {})'.format(filt, ccd))
            target_data_orig, target_mask = log.openData(ccd, '1', save=False,
                                                         mask=False)
            best_snr = 0
            best_snr_ap = '2'
            for ap in log.apnames[ccd]:

                if ap == '1':
                    continue
                comp_data, comp_mask = log.openData(ccd, ap, save=False,
                                                    mask=False)
                master_mask = np.logical_and(target_mask, comp_mask)
                target_data = target_data_orig[master_mask]
                comp_data = comp_data[master_mask]
                t_t, t_te, t_y, t_ye, _, _ = ((map(np.squeeze, np.split(target_data, 6, axis=1))))
                _, _, c_y, c_ye, _, _ = ((map(np.squeeze, np.split(comp_data, 6, axis=1))))
                diffFlux = t_y / c_y
                diffFluxErr = ((t_ye / t_y)**2 + (c_ye / c_y)**2)**0.5 * diffFlux
                comp_flux, comp_flux_err, snr, airmass = self.__cal_comp__(comp_data, filt, log.target_coords)
                calFlux = diffFlux * comp_flux.value
                calFluxErr = diffFluxErr * comp_flux.value
                # if ap == '2' and filt == 'u':
                #     start, end = self.__get_eclipse__(t_t, diffFlux, width=1)

                _, comp_mag_err = flux_to_ABmag(comp_flux, comp_flux_err)
                print('Aperture {} SNR = {:.2f}, Flux cal err = {:.3f} '
                      'mags = {:.3f}% (airmass = {:.2f})'.format(ap, snr, comp_mag_err,
                                              comp_flux_err*100/comp_flux, airmass))
                # print('{} band, aperture {}, flux calibration uncertainty = {:.3f} mags = {:.3f}%'.format(filt, ap, comp_mag_err, comp_flux_err*100/comp_flux))
                t_out = t_t[(t_t > start) & (t_t < end)]
                exp_out = t_te[(t_t > start) & (t_t < end)]
                calFlux_out = calFlux[(t_t > start) & (t_t < end)]
                calFluxErr_out = calFluxErr[(t_t > start) & (t_t < end)]
                weights = np.ones(len(t_out))
                bmjd_tdb = t_out + log.barycorr(t_out).value

                _, ax = plt.subplots()
                ax.scatter(bmjd_tdb, calFlux_out)
                plt.show()

                out = np.column_stack((bmjd_tdb, exp_out, calFlux_out,
                                       calFluxErr_out, weights, weights))
                fname = '{}_{}_{}_ap{}_fc.dat'.format(target, log.run, filt, ap)
                fname = os.path.join(log.path, 'reduced', fname)
                out_dict['data'][ap] = out
                out_dict['fname'][ap] = fname

                if snr > best_snr:
                    best_snr = snr
                    best_snr_ap = ap
                
            save_ap = input('Aperture to save [{}]: '.format(best_snr_ap))
            if not save_ap:
                save_ap = best_snr_ap
            if not os.path.isdir(os.path.join(os.path.join(log.path, 'reduced'))):
                os.makedirs(os.path.join(os.path.join(log.path, 'reduced')))
            np.savetxt(out_dict['fname'][save_ap], out_dict['data'][save_ap],
                       fmt='%11.9f %9.4e %9.4e %9.4e %1.0f %1.0f')


if __name__ == "__main__":
    # obs = Observation('hipercam')
    obs = Observation('ultracam')
    
    GD108_mags = dict(mean=dict(u=13.198, g=13.328, r=13.772, i=14.099, z=14.431),
                       err=dict(u=0.02, g=0.02, r=0.02, i=0.02, z=0.02))
    SA100_280_mags = dict(mean=dict(u=13.147-0.06, g=12.000+0.0049, r=11.691-0.003, i=11.606-0.0003, z=11.605+0.0032),
                          err=dict(u=0.02, g=0.02, r=0.02, i=0.02, z=0.02))
    Feige_67_mags = dict(mean=dict(u=11.0252, g=11.5116, r=12.0786, i=12.4567, z=13.4091),
                          err=dict(u=0.02, g=0.02, r=0.02, i=0.02, z=0.02))
    EG131_mags = dict(mean=dict(u=12.228, g=12.2196, r=12.359, i=12.544, z=12.750),
                          err=dict(u=0.02, g=0.02, r=0.02, i=0.02, z=0.02))
    WD2309_mags = dict(mean=dict(u=12.395, g=12.812, r=13.323, i=13.733, z=14.079),
                          err=dict(u=0.02, g=0.02, r=0.02, i=0.02, z=0.02))
    Feige110_mags = dict(mean=dict(u=11.168, g=11.547, r=12.059, i=12.446, z=12.787),
                          err=dict(u=0.02, g=0.02, r=0.02, i=0.02, z=0.02))
    GD71_mags = dict(mean=dict(u=12.451, g=12.775, r=13.277, i=13.651, z=14.034),
                     err=dict(u=0.02, g=0.02, r=0.02, i=0.02, z=0.02))
    GD153_mags = dict(mean=dict(u=12.692, g=13.072, r=13.594, i=13.984, z=14.361),
                     err=dict(u=0.02, g=0.02, r=0.02, i=0.02, z=0.02))
    LTT2415_mags = dict(mean=dict(u=13.156, g=12.354, r=12.110, i=12.033, z=12.053),
                     err=dict(u=0.02, g=0.02, r=0.02, i=0.02, z=0.02))
    EG274_mags = dict(mean=dict(u=10.699, g=10.804, r=11.255, i=11.580, z=11.922),
                     err=dict(u=0.02, g=0.02, r=0.02, i=0.02, z=0.02))
    
    # obs.add_observation(name='PPMXL', logfiles=['ultracam/2022_03_04/run026.log'], obs_type='atm')
    # obs.add_observation(name='CRTSJ0357', logfiles=['ultracam/2022_03_04/run020.log'], obs_type='atm')
    # obs.get_atm_ex()
    # obs.add_observation(name='GD71', logfiles=['ultracam/2022_03_04/run017.log'], obs_type='std', cal_mags=GD71_mags)
    # obs.get_zeropoint()
    # obs.add_observation(name='LTT2415', logfiles=['ultracam/2022_03_04/run024.log'], obs_type='std', cal_mags=LTT2415_mags)
    # obs.get_zeropoint()
    # obs.add_observation(name='EG274', logfiles=['ultracam/2022_03_04/run033.log'], obs_type='std', cal_mags=EG274_mags)
    # obs.get_zeropoint()
    # obs.add_observation(name='SDSSJ0624', logfiles=['ultracam/2022_03_04/run021.log'], obs_type='science')
    # obs.calibrate_science('SDSSJ0624', eclipse=False)

    # obs.clear()

    # obs.add_observation(name='1716bx', logfiles=['ultracam/2022_03_05/run034.log'], obs_type='atm')
    # # obs.add_observation(name='CRTSJ0357', logfiles=['ultracam/2022_03_05/run020.log'], obs_type='atm')
    # obs.get_atm_ex()
    # obs.add_observation(name='ZTFJ0752', logfiles=['ultracam/2022_03_05/run022.log'], obs_type='atm')
    # obs.get_atm_ex()
    # obs.add_observation(name='GD71', logfiles=['ultracam/2022_03_05/run017.log'], obs_type='std', cal_mags=GD71_mags)
    # obs.get_zeropoint()
    # obs.add_observation(name='EG274', logfiles=['ultracam/2022_03_05/run037.log'], obs_type='std', cal_mags=EG274_mags)
    # obs.get_zeropoint()
    # obs.add_observation(name='SDSSJ0624', logfiles=['ultracam/2022_03_05/run019.log'], obs_type='science')
    # obs.calibrate_science('SDSSJ0624', eclipse=False)

    # obs.clear()

    obs.add_observation(name='SDSSJ0624_atm', logfiles=['ultracam/2022_03_07/run030_atm.log'], obs_type='atm')
    # obs.add_observation(name='CRTSJ0357', logfiles=['ultracam/2022_03_07/run020.log'], obs_type='atm')
    obs.get_atm_ex()
    obs.add_observation(name='ZTFJ1802', logfiles=['ultracam/2022_03_07/run040.log'], obs_type='atm')
    obs.get_atm_ex()
    obs.add_observation(name='GD71', logfiles=['ultracam/2022_03_07/run019.log'], obs_type='std', cal_mags=GD71_mags)
    obs.get_zeropoint()
    obs.add_observation(name='GD108', logfiles=['ultracam/2022_03_07/run024.log'], obs_type='std', cal_mags=GD108_mags)
    obs.get_zeropoint()
    obs.add_observation(name='GD153', logfiles=['ultracam/2022_03_07/run033.log'], obs_type='std', cal_mags=GD153_mags)
    obs.get_zeropoint()
    obs.add_observation(name='SDSSJ0624', logfiles=['ultracam/2022_03_07/run030.log'], obs_type='science')
    obs.calibrate_science('SDSSJ0624', eclipse=False)