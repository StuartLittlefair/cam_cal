import numpy as np
import numpy.ma as ma
import astropy.units as u


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


def mask_data(t_data, t_mask, c_data, c_mask):
    negative_mask = ma.getmask(ma.masked_greater(c_data[:,2], 0))
    master_mask = np.logical_and.reduce((t_mask, c_mask, negative_mask))
    target_data = t_data[master_mask]
    comp_data = c_data[master_mask]
    return target_data, comp_data


def top_tail(arrays, index_array, lolim, uplim):
        array = np.column_stack(arrays)
        cut_array = array[(index_array > lolim) & (index_array < uplim)]
        return [array.flatten() for array in np.split(cut_array, cut_array.shape[1], axis=1)]


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