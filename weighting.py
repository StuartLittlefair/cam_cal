import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, least_squares, newton
from scipy.interpolate import interp1d
from scipy.integrate import simps
import os


def sigmoid_corr(x, x0, slope):
    """
    A sigmoid function with a correction term so that x0 reflects the
    start of the top of the sigmoid function rather than the centre of
    the slope.
    """
    exponent = slope * (-x + x0 - 4/slope)
    exponent[exponent > 700] = 700
    exponent[exponent < -700] = -700
    return 1 / (1 + np.exp(exponent))


def weight_boost(x, max, t_start, t_end, slope):
    weight_multplier = max * (sigmoid_corr(x, t_start, slope)
                        * sigmoid_corr(x, t_end, -slope))
    return weight_multplier


def in_eg_weights_boost(x, t1, t2, t3, t4, slope, scale):
    weights = np.ones_like(x)
    weights += weight_boost(x, scale, t1, t2, slope)
    weights += weight_boost(x, scale, t3, t4, slope)
    return weights


def weight_proportion(scale, slope, t, t1, t2, t3 ,t4):
    weights = in_eg_weights_boost(t, t1, t2, t3, t4, slope, scale)
    mask = ((t > t1) & (t < t2)) | ((t > t3) & (t < t4))
    area_whole = simps(weights)
    area_boosted = simps(weights[mask])
    if not mask.any():
        area_boosted = 0
    prop = area_boosted / area_whole
    return prop


def scale_weights(chosen_proportion, slope, t, t1, t2, t3 ,t4):
    fn = lambda scale ,slope, t, t1, t2, t3 ,t4, chosen_proportion: weight_proportion(scale, slope, t, t1, t2, t3 ,t4) - chosen_proportion
    res = least_squares(fn, x0=[5], args=[slope, t, t1, t2, t3 ,t4, chosen_proportion], method='trf', bounds=(0.1, 1000))
    return res.x


def get_weights(data, t1, t2, t3, t4, slope):
    t, exp, y, ye, _, esubd = data.T
    weightboost = in_eg_weights_boost(t, t1, t2, t3, t4, slope, 10)
    print('Scaling weights')
    scale = scale_weights(0.333, slope, t, t1, t2, t3, t4)
    print(f"Weights scaled by {float(scale):.2f}")
    print("Applying weights to light curves...")
    weight = in_eg_weights_boost(t, t1, t2, t3, t4, slope, scale)
    out = np.column_stack([t, exp, y, ye, weight, esubd])
    return out
