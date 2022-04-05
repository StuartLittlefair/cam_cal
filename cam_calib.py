import numpy as np
import pandas as pd
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

from logfile import Logfile
import utils


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
        # set instrument specific variables
        if self.instrument=='ultracam':
            self.tel_location = EarthLocation.of_site('La Silla Observatory')
            self.filt2ccd = {'u':'3', 'us':'3', 'g':'2', 'gs':'2', 'r':'1',
                             'rs':'1', 'i':'1', 'is':'1', 'z':'1', 'zs':'1'}
            self.rootDataDir = '/local/alex/backed_up_on_astro3/Data/photometry/ultracam'
            self.stds = pd.read_csv('cam_standards/ucam_flux_stds.csv', dtype=format_dict)
        elif self.instrument=='hipercam':
            self.tel_location = EarthLocation.of_site('Roque de los Muchachos')
            self.filt2ccd = {'us':'1', 'gs':'2', 'rs':'3', 'is':'4', 'zs':'5'}
            self.rootDataDir = '/local/alex/backed_up_on_astro3/Data/photometry/hipercam'
            self.stds = pd.read_csv('cam_standards/hcam_flux_stds.csv', dtype=format_dict)
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
                                  airmass={'airmass':0})
        elif self.instrument=='hipercam':
            self.zeropoint = dict(mean={'us':28.15, 'gs':29.22, 'rs':28.78, 'is':28.43, 'zs':27.94},
                                  err={'us':0.05, 'gs':0.05, 'rs':0.05, 'is':0.05, 'zs':0.05},
                                  airmass={'airmass':0})


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


    def __atm_chisq__(self, pars, apertures):
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
        # TODO: Either prevent looping through names or average zeropoints.
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
        mag_i, mag_i_err = utils.flux_to_mag(flux, flux_err)
        mag_i0 = mag_i - (self.atm_extinction['mean'][filt] * airmass)
        mag_comp = mag_i0 + self.zeropoint['mean'][filt]
        mag_comp_err = self.__flux_cal_err__(airmass, mag_i_err, filt)
        snr_comp = np.median(flux/flux_err)
        comp_flux, comp_flux_err = utils.magAB_to_flux(np.mean(mag_comp), np.mean(mag_comp_err))
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
    

    def calibrate_science(self, target_name, eclipse=None):
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
        negative_mask = ma.getmask(ma.masked_greater(comp_data[:,2], 0))
        master_mask = np.logical_and.reduce((target_mask, comp_mask, negative_mask))
        target_data = target_data_orig[master_mask]
        comp_data = comp_data[master_mask]
        t_t, t_te, t_y, t_ye, _, _ = ((map(np.squeeze, np.split(target_data, 6, axis=1))))
        _, _, c_y, c_ye, _, _ = ((map(np.squeeze, np.split(comp_data, 6, axis=1))))
        diffFlux = t_y / c_y

        if eclipse:
            start, end = self.__get_eclipse__(t_t, diffFlux, width=eclipse)
        else:
            start, end = t_t[0]-0.0001, t_t[-1]+0.0001 

        for filt in log.filters:
            out_dict = dict(data=dict(), fname=dict())
            ccd = self.filt2ccd[filt]
            print(f"\n{filt}-band (CCD {ccd})")
            target_data_orig, target_mask = log.openData(ccd, '1', save=False,
                                                         mask=False)
            best_snr = 0
            best_snr_ap = '2'
            for ap in log.apnames[ccd]:

                if ap == '1':
                    continue
                comp_data, comp_mask = log.openData(ccd, ap, save=False,
                                                    mask=False)
                negative_mask = ma.getmask(ma.masked_greater(comp_data[:,2], 0))
                master_mask = np.logical_and.reduce((target_mask, comp_mask, negative_mask))
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

                _, comp_mag_err = utils.flux_to_ABmag(comp_flux, comp_flux_err)
                comp_err_percent = comp_flux_err*100/comp_flux
                print(f"Aperture {ap} SNR = {snr:.2f}, Flux cal err = {comp_mag_err:.3f} "
                      f"mags = {comp_err_percent:.3f}% (airmass = {airmass:.2f})")

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
                fname = f"{target}_{log.run}_{filt}_ap{ap}_fc.dat"
                fname = os.path.join(log.path, 'reduced', fname)
                out_dict['data'][ap] = out
                out_dict['fname'][ap] = fname
                

                if snr > best_snr:
                    best_snr = snr
                    best_snr_ap = ap
                
            save_ap = input(f"Aperture to save [{best_snr_ap}]: ")
            if not save_ap:
                save_ap = best_snr_ap
            if not os.path.isdir(os.path.join(os.path.join(log.path, 'reduced'))):
                os.makedirs(os.path.join(os.path.join(log.path, 'reduced')))
            np.savetxt(out_dict['fname'][save_ap], out_dict['data'][save_ap],
                       fmt='%11.9f %9.4e %9.4e %9.4e %1.0f %1.0f')


if __name__ == "__main__":
    # obs = Observation('hipercam')
    obs = Observation('ultracam')
    
    Feige_67_mags = dict(mean=dict(u=11.0252, g=11.5116, r=12.0786, i=12.4567, z=13.4091),
                          err=dict(u=0.02, g=0.02, r=0.02, i=0.02, z=0.02))

    obs.add_observation(name='CXOUJ110926.4-650224', logfiles=['/local/alex/backed_up_on_astro3/Data/photometry/ultracam/2022_03_02/run010.log'], obs_type='atm')
    obs.get_atm_ex()
    obs.add_observation(name='GD153', logfiles=['/local/alex/backed_up_on_astro3/Data/photometry/ultracam/2022_03_02/run011.log'], obs_type='std', cal_mags='GD153')
    obs.get_zeropoint()
    obs.add_observation(name='1712af', logfiles=['/local/alex/backed_up_on_astro3/Data/photometry/ultracam/2022_03_02/run015.log'], obs_type='science')
    obs.calibrate_science('1712af', eclipse=1.5)