import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import heapq


def sigmoid(x, x0, slope):
    """
    Sigmoid function with limited exponents.
    """
    exponent = slope * (-x + x0)
    exponent[exponent > 700] = 700
    exponent[exponent < -700] = -700
    return 1 / (1 + np.exp(exponent))


def sigmoid_fit_func(x, scale, minimum, slope, offset1, offset2, A, t0):
    """
    A quick and dirty eclipse-like function made up of two opposing sigmoids
    and a quadratic term centred on mid-eclipse.
    """
    func = scale * (sigmoid(x, offset1, -slope) + sigmoid(x, offset2, slope)) + minimum + A*(x - t0)**2
    return func


def sigmoid_residuals(pars, t, y):
    return y - sigmoid_fit_func(t, *pars)


def get_ingress_egress(time, flux):
    """Roughly estimate times of ingress and egress."""
 
    # half_flux = (np.max(flux) + np.min(flux)) / 2
    half_flux = (np.median(heapq.nlargest(10, flux)) + np.median(heapq.nsmallest(10, flux))) / 2
    a1 = flux.copy()
    idx1 = np.nanargmin(np.abs(a1 - half_flux))
    a1[idx1-4:idx1+5] = np.nan
    idx2 = np.nanargmin(np.abs(a1 - half_flux))
    interp1 = interp1d(flux[idx1-4:idx1+5], time[idx1-4:idx1+5])
    interp2 = interp1d(flux[idx2-4:idx2+5], time[idx2-4:idx2+5])
    t1 = interp1(half_flux)
    t2 = interp2(half_flux)
    return t1, t2


def contact_points(t, y):
    """
    Estimate eclipse contact points by fitting a quick and dirty eclipse-like
    function made up of two opposing sigmoids and a quadratic term centred on
    mid-eclipse (to account for reflection effect/ellipsoidal modulation).
    """
    smoothed = savgol_filter(y, 20, 3)
    t1, t2 = get_ingress_egress(t, smoothed)
    t0 = np.mean([t1, t2])
    minimum = np.abs(np.median(heapq.nsmallest(50, smoothed)))
    top = heapq.nlargest(len(y)-50, smoothed)
    depth = np.median(top) - minimum
    res = least_squares(sigmoid_residuals, x0=[depth, minimum, 10000, t1, t2, 0, t0], method='lm', args=[t, y])
    pars = res.x.copy()
    res.x[-2] = 0
    tdense = np.linspace(t[0], t[-1], 1000)
    fit_depth = np.max(sigmoid_fit_func(tdense, *res.x)) - np.min(sigmoid_fit_func(tdense, *res.x))
    t_half_1 = tdense[tdense < t0]
    t_half_2 = tdense[tdense > t0]
    rev_interp_1 = interp1d(sigmoid_fit_func(t_half_1, *res.x), t_half_1)
    rev_interp_2 = interp1d(sigmoid_fit_func(t_half_2, *res.x), t_half_2)
    t1 = rev_interp_1(0.98*fit_depth + np.min(sigmoid_fit_func(tdense, *res.x)))
    t2 = rev_interp_1(0.02*fit_depth + np.min(sigmoid_fit_func(tdense, *res.x)))
    t3 = rev_interp_2(0.02*fit_depth + np.min(sigmoid_fit_func(tdense, *res.x)))
    t4 = rev_interp_2(0.98*fit_depth + np.min(sigmoid_fit_func(tdense, *res.x)))
    return t1, t2, t3, t4