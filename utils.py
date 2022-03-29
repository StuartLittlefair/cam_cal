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