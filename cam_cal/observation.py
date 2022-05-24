import numpy as np
import pandas as pd
from hipercam import hlog
import matplotlib.pyplot as plt
from astropy.table import Table, QTable
import numpy.ma as ma
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, name_resolve
import astropy.units as u
from astropy.time import Time
import time
from tkinter.filedialog import askopenfilenames
from tkinter import Tk
import os
import re
import warnings
from mergedeep import merge
from pkg_resources import resource_filename

from cam_cal.fits import FITS
from cam_cal.logfile import Logfile
import cam_cal.utils as utils
import cam_cal.times as times
import cam_cal.weighting as weighting

fpath = resource_filename('cam_cal', 'cam_standards/')

warnings.filterwarnings('error')


class Observation:
    """Class for reducing and flux calibrating ULTRACAM and HiPERCAM
       eclipse light curves."""


    def __init__(self, instrument='ultracam', tel_location=None):
        root = Tk()
        root.withdraw()
        self.instrument = instrument
        format_dict = {'Name': 'string',
                       'RA_hms': 'string',
                       'DEC_dms': 'string',
                       'Type': 'string',
                       'SpType': 'string',
                       'fname': 'string',
                       'Variability': 'string'}
        self.header_dict = dict()
        # set instrument specific variables
        if self.instrument=='ultracam':
            self.tel_location = EarthLocation.of_site('La Silla Observatory')
            self.filt2ccd = {'u':'3', 'us':'3', 'g':'2', 'gs':'2', 'r':'1',
                             'rs':'1', 'i':'1', 'is':'1', 'z':'1', 'zs':'1'}
            self.rootDataDir = os.environ.get("UCAM_DATA", "/home")
            self.stds = pd.read_csv(f"{fpath}ucam_flux_stds.csv", dtype=format_dict)
        elif self.instrument=='hipercam':
            self.tel_location = EarthLocation.of_site('Roque de los Muchachos')
            self.filt2ccd = {'us':'1', 'gs':'2', 'rs':'3', 'is':'4', 'zs':'5'}
            self.rootDataDir = os.environ.get("HCAM_DATA", "/home")
            self.stds = pd.read_csv(f"{fpath}hcam_flux_stds.csv", dtype=format_dict)
        else: raise ValueError(f"{self.instrument} is not a valid instrument")
        if tel_location:
            if tel_location in EarthLocation.get_site_names():
                self.tel_location = EarthLocation.of_site(tel_location)
            elif isinstance(tel_location, EarthLocation):
                self.tel_location = tel_location
            else:
                raise ValueError(f'"{tel_location}" not a valid EarthLocation site name or '
                                 'an EarthLocation object')
        # set extinction coeffs and instrumental zeropoints
        self.set_default()

        #initialise observation dict
        self.observations = dict.fromkeys(['science', 'std', 'atm'])

    
    def set_default(self):
        """Sets default atmospheric extinctions and instrument zeropoints"""

        # self.atm_extinction = dict(mean=dict(u=0.48, g=0.24, r=0.18, i=0.15, z=0.10),
        #                            err=dict(u=0.05, g=0.05, r=0.05, i=0.05, z=0.05))
        self.atm_extinction = dict(mean={'us':0.48, 'gs':0.24, 'rs':0.18, 'is':0.15, 'zs':0.10},
                                   err={'us':0.05, 'gs':0.05, 'rs':0.05, 'is':0.05, 'zs':0.05})
        if self.instrument=='ultracam':
            self.zeropoint = dict(mean={'us':25.09, 'gs':26.70, 'rs':26.30, 'is':25.90, 'zs':25.30},
                                  err={'us':0.05, 'gs':0.05, 'rs':0.05, 'is':0.05, 'zs':0.05},
                                  airmass=0.0)
        elif self.instrument=='hipercam':
            self.zeropoint = dict(mean={'us':28.15, 'gs':29.22, 'rs':28.78, 'is':28.43, 'zs':27.94},
                                  err={'us':0.05, 'gs':0.05, 'rs':0.05, 'is':0.05, 'zs':0.05},
                                  airmass=0.0)


    def check_path(self, fname):
        import errno
        if not os.path.exists(fname):
            fname = os.path.join(self.rootDataDir, fname)
        if not os.path.exists(fname):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), fname)
        return fname
            


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
        else:
            logfiles = [self.check_path(logfile) for logfile in logfiles]

        if obs_type=='std' and not cal_mags:
            raise ValueError('If observations are of a flux standard then '
                             'calibrated magnitudes must be supplied')
        obs1, obs2, obs3 = dict(), dict(), dict()
        obs1['name'], obs1['logfiles'] = name, logfiles
        if obs_type == 'std':
            obs1['cal_mags'] = self.match_std(cal_mags)
        obs2[name] = obs1
        obs3[obs_type] = obs2
        self.observations = merge(self.observations, obs3)

    
    def construct_dict(self, series):
        mags_dict = dict(mean={'us':0, 'gs':0, 'rs':0, 'is':0, 'zs':0},
                         err={'us':0.02, 'gs':0.02, 'rs':0.02, 'is':0.02, 'zs':0.02})
        for filt in ['us', 'gs', 'rs', 'is', 'zs']:
            mags_dict['mean'][filt] = float(series[filt])
        return mags_dict


    def match_std(self, standard):
        if isinstance(standard, str):
            if self.stds['Name'].isin([standard]).any():

                return self.construct_dict(self.stds.loc[self.stds['Name'] == standard])
            else:
                try:
                    standard = SkyCoord.from_name(standard, parse=True)
                except name_resolve.NameResolveError:
                    raise ValueError(f"{standard} doesn't match the list of standards and can't be resolved in SIMBAD.")

        if isinstance(standard, SkyCoord):
            catalog = SkyCoord(self.stds['RA']*u.deg, self.stds['DEC']*u.deg)
            idx, d2d, d3d = standard.match_to_catalog_sky(catalog)
            return self.construct_dict(self.stds.iloc[idx])

        elif isinstance(standard, dict):
            return standard

        else: raise ValueError(f"'{standard}' is not a valid input.")
    

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
        #TODO: add theilslopes method
        
        if 'atm' not in self.observations.keys():
            raise ValueError('No atmospheric extinction observations added.')
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
            p, errs = self.fit_atm_ext(filt, apertures)
            print(f"{filt}-band extinction: {p[0]:.3f} +- {errs[0]:.3f}")
            if plot:
                _, ax = plt.subplots()
                for idx, apdata in enumerate(apertures):
                    ax.scatter(apdata[:,0], utils.flux_to_mag(apdata[:,2], apdata[:,3])[0],
                               label="ap{}, {}".format(*ap_info[idx]))
                    ax.plot(apdata[:,0], utils.linear(apdata[:,0], p[0], p[idx+1]), 'k--')
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


    def atm_chisq(self, pars, apertures):
        """Computes combined chi (not squared) for every aperture"""

        k, *c = pars
        diff = np.array([])
        for i, data in enumerate(apertures):
            # compute chi residuals for each aperture
            mag, mag_err = utils.flux_to_mag(data[:,2], data[:,3])
            airmass = data[:,0]
            diff_new = np.abs((mag - utils.linear(airmass, k, c[i])) / mag_err)
            diff = np.hstack((diff, diff_new))
        return diff
    

    def red_chisq(self, pars, apertures):
        """Reduced chi-squared given data and fit parameters"""

        chisq = np.sum(self.atm_chisq(pars, apertures)**2)
        ndata = 0
        for data in apertures:
            ndata += data.shape[0]
        ndf = ndata-len(pars)
        red_chisq = chisq / ndf
        return red_chisq.value


    def fit_atm_ext(self, filt, apertures):
        """Fits atmospheric extinction coefficient to data in given apertures."""
        
        # add offset parameter for each aperture
        x0 = [self.atm_extinction['mean'][filt]] + [0 for i in apertures]
        p, cov, _, _, _ = leastsq(self.atm_chisq, x0,
                                  full_output=True, args=(apertures))
        red_chisq = self.red_chisq(p, apertures)
        # scale variances with reduced chisqr and get std_err for fitted params
        errs = np.sqrt(np.diag(cov * red_chisq))
        return p, errs


    def get_zeropoint(self):
        """Calculate magnitude of star that would give 1 count/sec at the
           instrument if no atmosphere. Must have added a 'std' observation
           with calibrated magnitudes in the filter system used."""
        # TODO: Either prevent looping through names or average zeropoints.
        # Need to check multiple logfiles for same target at different airmasses
        # check zip(target_names, logfiles) will always be the same length

        if 'std' not in self.observations.keys():
            raise ValueError('No flux standard observations added.')
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
                mag_i, mag_i_err = utils.flux_to_mag(flux, flux_std_err)
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
                zp_dict['airmass'] = np.mean(airmass.value)
                print(f"{filt}-band zeropoint: {zp:.3f} +- {zp_err:.3f} (airmass = {np.mean(airmass):.2f})")
        write = input("Set calculated zeropoints [y/n]? ")
        if not write or write=='y':
            for key, value in zp_dict.items():
                self.zeropoint[key] = value
            self.standard = log.target
            self.std_run = log.run
            print('Zeropoints set.\n')
        else:
            print('Keeping default zeropoints.\n')

    
    def flux_cal_err(self, airmass_t, mag_i_err, filt):
        mstd_mstdi_var = (self.zeropoint['err'][filt]**2
                          - (self.zeropoint['airmass']
                             * self.atm_extinction['err'][filt])**2)
        am_diff = self.zeropoint['airmass'] - airmass_t
        airmass_var = (self.atm_extinction['err'][filt] * am_diff)**2
        mag_t_err = np.sqrt( mstd_mstdi_var + airmass_var +  mag_i_err**2 )
        return mag_t_err


    def cal_comp(self, data, filt, coords):
        """Calibrate the mean flux of the comparison and it's uncertainty using
           the instrumental zeropoint and atmospheric extinction"""
        # TODO: implement outputting/storing comparison magnitude.

        t, te, y, ye, _, _ = data.T
        flux = y/(te*86400)
        flux_err = ye/(te*86400)
        t = Time(t, format='mjd', scale='tdb')
        airmass = self.airmass(t, coords)

        mag_i, mag_i_err = utils.flux_to_mag(flux, flux_err)
        mag_i0 = mag_i - (self.atm_extinction['mean'][filt] * airmass)
        mag_comp = mag_i0 + self.zeropoint['mean'][filt]
        mag_comp_err = self.flux_cal_err(airmass, mag_i_err, filt)
        comp_flux, comp_flux_err = utils.magAB_to_flux(np.mean(mag_comp), np.mean(mag_comp_err))
        return comp_flux, comp_flux_err, np.mean(airmass)


    def get_eclipse(self, t1, t4, width=1):
        """Give times required to chop out width x eclipsewidth either side of the middle of eclipse.
           width=1 gives the eclipse plus half the eclipse width either side."""

        eclipsewidth = t4 - t1
        t0 = (t1 + t4) / 2
        start = t0 - (width * eclipsewidth)
        end = t0 + (width * eclipsewidth)
        return start, end


    def diff_phot(self, log, ccd, ap, target_data, target_mask):
        comp_data, comp_mask = log.openData(ccd, ap, save=False,
                                            mask=False)
        target, comp = utils.mask_data(target_data, target_mask, comp_data, comp_mask)

        _, _, t_y, t_ye, _, _ = target.T
        _, _, c_y, c_ye, _, _ = comp.T
        comp_snr = np.median(c_y/c_ye)
        diffFlux = t_y / c_y
        diffFluxErr = ((t_ye / t_y)**2 + (c_ye / c_y)**2)**0.5 * np.abs(diffFlux)
        return target, comp, diffFlux, diffFluxErr, comp_snr


    def calibrate_science(self, target_name, comp_mag=None, comp_mag_err=None, eclipse=None):
        """
        Flux calibrate the selected science target using the calibrated comparison stars.
        eclipse width can be specified.

        """

        data_arrays = dict()
        log = Logfile(self.observations['science'][target_name]['logfiles'][0],
                      self.instrument, self.tel_location)
        target = log.target.replace(' ', '_')
        ccd = self.filt2ccd['g']
        target_data_orig, target_mask = log.openData(ccd, '1', save=False,
                                                     mask=False)
        target_data, _, diffFlux, _, _ = self.diff_phot(log, ccd, '2', target_data_orig, target_mask)
        t_t, t_te, _, _, _, _ = target_data.T

        if eclipse:
            t1, t2, t3, t4 = times.contact_points(t_t, diffFlux)
            start, end = self.get_eclipse(t1, t4, width=eclipse)

        else:
            start, end = t_t[0]-0.0001, t_t[-1]+0.0001

        data_arrays = dict(data=dict(), header=dict())
        for filt in log.filters:
            out_dict = dict()
            ccd = self.filt2ccd[filt]
            print(f"\n{filt}-band (CCD {ccd})")
            target_data_orig, target_mask = log.openData(ccd, '1', save=False,
                                                         mask=False)
            best_snr = 0
            best_snr_ap = '2'
            for ap in log.apnames[ccd]:

                if ap == '1':
                    continue

                target_data, comp_data, diffFlux, diffFluxErr, comp_snr = self.diff_phot(log, ccd, ap, target_data_orig, target_mask)
                t_t, t_te, _, _, _, _ = target_data.T

                
                if comp_mag and comp_mag_err:
                    comp_flux, comp_flux_err = utils.magAB_to_flux(comp_mag, comp_mag_err)
                    airmass = 0.00
                else:
                    comp_flux, comp_flux_err, airmass = self.cal_comp(comp_data, filt, log.target_coords)

                calFlux = diffFlux * comp_flux.value
                calFluxErr = diffFluxErr * comp_flux.value

                _, comp_mag_err = utils.flux_to_ABmag(comp_flux, comp_flux_err)
                comp_err_percent = comp_flux_err*100/comp_flux
                print(f"Aperture {ap} SNR = {comp_snr:.2f}, Flux cal err = {comp_mag_err:.3f} "
                      f"mags = {comp_err_percent:.3f}% (airmass = {airmass:.2f})")

                t_out, exp_out, calFlux_out, calFluxErr_out = utils.top_tail([t_t, t_te, calFlux, calFluxErr], t_t, start, end)
                weights = np.ones(len(t_out))
                bary_corr = log.barycorr(t_out).value
                bmjd_tdb = t_out + bary_corr

                _, ax = plt.subplots()
                ax.scatter(bmjd_tdb, calFlux_out)
                plt.show()

                out = np.column_stack((bmjd_tdb, exp_out, calFlux_out,
                                       calFluxErr_out, weights, weights))
                if eclipse:
                    slope = 150000 / np.median(np.diff(bmjd_tdb * 86400))
                    out = weighting.get_weights(out, t1+bary_corr, t2+bary_corr, t3+bary_corr, t4+bary_corr, slope)
                out_dict[ap] = out


                if comp_snr > best_snr:
                    best_snr = comp_snr
                    best_snr_ap = ap

            save_ap = input(f"Aperture to save [{best_snr_ap}]: ")
            if not save_ap:
                save_ap = best_snr_ap
            hdr = dict(FILTER=filt,
                       COMP_AP=save_ap,
                       ATM_EX=self.atm_extinction['mean'][filt],
                       ATM_EX_E=self.atm_extinction['err'][filt],
                       ZP=self.zeropoint['mean'][filt],
                       ZP_E=self.zeropoint['err'][filt])

            if not os.path.isdir(os.path.join(log.path, 'reduced')):
                os.makedirs(os.path.join(log.path, 'reduced'))
            if not os.path.isdir(os.path.join(log.path, 'reduced', target)):
                os.makedirs(os.path.join(log.path, 'reduced', target))

            tab_data = Table(out_dict[save_ap], names=('BMJD(TDB)', 'exp_time', 'flux', 'flux_err', 'weight', 'esubd'))
            data_arrays['data'][filt] = tab_data
            data_arrays['header'][filt] = hdr
        t = Time(start, format='mjd', scale='tdb').ymdhms
        if t[3] < 14:
            night = f"{t[0]}-{t[1]}-{t[2]-1}"

        header = dict(TARGET=target,
                      RUN=log.run, FILTERS=", ".join(log.filters),
                      NIGHT=night,
                      INSTR=self.instrument,
                      RA=log.target_coords.ra.deg,
                      DEC=log.target_coords.dec.deg,
                      STD_STAR=self.standard,
                      STD_RUN=self.std_run,
                      STD_AIRM=self.zeropoint['airmass'],
                      ZP=", ".join([str(val) for val in list(self.zeropoint['mean'].values())]),
                      ZP_E=", ".join([str(val) for val in list(self.zeropoint['err'].values())]),
                      ATM_EX=", ".join([str(val) for val in list(self.atm_extinction['mean'].values())]),
                      ATM_EX_E=", ".join([str(val) for val in list(self.atm_extinction['err'].values())]),
                      TIME_CAL=time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
                      TIME_ID=int(time.time()*1000)
                      )
        
        hdul = FITS.create(data_arrays, header)
        fname = f"{target}.fits"
        path = os.path.join(log.path, 'reduced', target)
        fname = os.path.join(path, fname)
        hdul.write(fname)
        hdul.to_lcurve(filepath=path)
        # write_FITS(fname, data_arrays, header)