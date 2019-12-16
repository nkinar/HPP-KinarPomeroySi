import numpy as np
from numpy import pi
from ei import ei_inv_neg_array
from deriv_forward_inv import forward_first_derivative, inverse_first_derivative
from get_size import length
from ei import ei
from scipy.optimize import minimize
from forward_model_nominal import compute_delta_T_dual_heating_cooling, compute_delta_T_dual_infinite_line_heating
np.seterr(invalid='ignore')

"""
Here is the code to compute the signal processing dual probe (DP) model.
"""


##############################################################################################
from get_eps import get_eps
eps = get_eps()


def ei_inv_cooling(fxy, th, tc):
    """
    Function to compute the inverse of f(x, y) = Ei(-x) - Ei(-y)
    NOTE that Nelder-Mead function works well, and that an operation can cause underflow.
    :param fxy:     as the function
    :param th:      as the time of heating
    :param tc:      as the cooling time vector
    :return:
    """
    def _func(x, *args):
        hh = x[0]
        known = args[0]
        th_known = args[1]
        t_known = args[2]
        upper = args[3]
        xx = hh / (t_known-th_known)
        yy = hh / t_known
        first = ei(-xx)
        second = ei(-yy)
        calc = first - second   # this operation can be subject to invalid floating point math
        if calc < 0:
            calc = eps
        elif calc > upper:
            calc = upper-eps
        out = known-calc
        if out < 0:             # compute absolute value (do this directly in case multiple precision math is required)
            out = -out
        return out
    nn = length(fxy)
    recon = np.zeros(nn)
    hstart = 80
    xstart = (hstart,)
    th_known = th
    for k in range(nn):
        t_known = tc[k]
        term = t_known-th_known
        known = fxy[k]
        z = np.sqrt(term/t_known)
        upper = term*(((1+z)/z)*(1/known))**2
        args = (known, th, t_known, upper)
        res = minimize(_func, xstart, args, method='Nelder-Mead', tol=1.0e-30)  # NOTE THE TOLERANCE REQUIRED
        recon[k] = res.x[0]
    return recon


def dp_model_recomp_shared(k, t, q, gamma4, dt, r0):
    """
    Shared code to recompose the model
    :param k:                   as the thermal conductivity
    :param t:                   as the time vector
    :param q:                   as the heat input into the soil
    :param gamma4:              as the gamma4 processed outside of this function
    :param dt:                  as the timestep
    :param r0:                  as the known r0 (or r0=None if the r0 is to be computed)
    :return:
    """
    gamma5 = np.sqrt(np.abs(gamma4))
    gamma6 = np.log(gamma5)
    gamma7 = forward_first_derivative(gamma6, dt)
    log_rinv = np.log(r0)
    g8 = inverse_first_derivative(gamma7, dt, log_rinv)
    g9 = np.exp(g8)
    r_t = g9
    return r_t


##############################################################################################


def dp_model_recomp_heating_cooling_trim(fs, q, k, t, th, dT_known, r0, get_gamma4=False):
    """
    Run the DP model, but take into consideration a trimmed curve that can be:
    (1) heating and cooling
    (2) cooling only
    DO NOT use this function for heating only.  Use the dp_model_recomp_heating() function below if
    only the heating curve is passed in.
    :param fs:                      as the sampling rate
    :param q:                       as the heat input into the soil
    :param k:                       as the thermal conductivity
    :param t:                       as the time vector
    :param th:                      as the threshold time
    :param dT_known:                as the known change in temperature
    :param r0:                      as the initial radius
    :param get_gamma4:              True to return gamma4 as well
    :return: r_t or (r_t, gamma4)
    """
    n = length(t)
    if length(dT_known) != n:
        raise ValueError('dp_model_recomp_heating_cooling: the length of t must be the same as dT')
    term = -(4 * pi * k) / q
    dt = 1 / fs
    if t[0] >= th:  # cooling only
        tvec_heat = []
        tvec_cool = t
        dT_known_heat = []
        dT_known_cool = dT_known
        gamma4_heat = []
    else:           # heating and cooling
        nheat = int(np.ceil(fs*(th-t[0])))
        tvec_heat = t[:nheat]
        tvec_cool = t[nheat:]
        dT_known_heat = dT_known[:nheat]
        dT_known_cool = dT_known[nheat:]
        gamma2_heat = term * dT_known_heat
        gamma3_heat = ei_inv_neg_array(gamma2_heat) * tvec_heat
        gamma4_heat = -gamma3_heat
    # COOLING
    dT_cool_strip = -term*dT_known_cool
    gamma4_cool = ei_inv_cooling(dT_cool_strip, th, tvec_cool)
    gamma4 = np.abs(np.concatenate((gamma4_heat, gamma4_cool)))
    r_t = dp_model_recomp_shared(k, t, q, gamma4, dt, r0)
    if get_gamma4:
        return r_t, gamma4
    return r_t


def dp_model_recomp_heating_cooling(fs, q, k, t, th, dT_known, r0, get_gamma4=False):
    """
    DP inverse model for heating and cooling.
    USE THIS FUNCTION ONLY FOR THE FULL CURVE THAT IS NOT TRIMMED.

    NOTE that the dT_known has to be sufficiently smooth for the derivative operation
    to not be contaminated by noise

    :param fs:                      as the sampling rate
    :param q:                       as the heat input into the soil
    :param k:                       as the thermal conductivity
    :param t:                       as the time vector
    :param th:                      as the time of heating
    :param dT_known:                as the change in temperature
    :param r0:                      as an assumed change in temperature
    :param get_gamma4:              True to obtain the gamma4
    :return:
    """
    n = length(t)
    if length(dT_known) != n:
        raise ValueError('dp_model_recomp_heating_cooling: the length of t must be the same as dT')
    term = -(4*pi*k)/q
    dt = 1 / fs
    nn = int(np.ceil(fs*th))
    tvec_heat = t[:nn]
    tvec_cool = t[nn:]
    dT_known_heat = dT_known[:nn]
    dT_known_cool = dT_known[nn:]
    # HEATING
    gamma2_heat = term*dT_known_heat
    gamma3_heat = ei_inv_neg_array(gamma2_heat) * tvec_heat
    gamma4_heat = -gamma3_heat
    # COOLING
    dT_cool_strip = -term*dT_known_cool
    gamma4_cool = ei_inv_cooling(dT_cool_strip, th, tvec_cool)
    gamma4 = np.concatenate((gamma4_heat, gamma4_cool))
    r_t = dp_model_recomp_shared(k, t, q, gamma4, dt, r0)
    if get_gamma4:
        return r_t, gamma4
    return r_t


def dp_model_recomp_heating(fs, q, k, t, dT_known, r0, get_gamma4=False):
    """
    DP inverse model for heating only.  NOTE that the dT_known has to be
    sufficiently smooth for the derivative operation to not be contaminated by noise
    
    :param fs:                      as the sampling frequency (Hz)
    :param q:                       as the heat input into the soil
    :param k:                       as the thermal conductivity
    :param t:                       as the time vector
    :param dT_known:                as the change in temperature
    :param r0:                      as the known r0
    :param get_gamma4:              True to return gamma4
    :return:
    """
    dt = 1 / fs
    term = -(4*pi*k)/q
    gamma2 = term*dT_known
    gamma3 = ei_inv_neg_array(gamma2) * t
    gamma4 = -gamma3
    r_t = dp_model_recomp_shared(k, t, q, gamma4, dt, r0)
    if get_gamma4:
        return r_t, gamma4
    return r_t


def dp_model_recomp_heating_unknown_r0(fs, q, k, alpha, t, dT_known, r0_start, get_gamma4=False):
    """
    DP inverse model for heating only when the r0 is not exactly known.  The r0_start
    is simply a starting value for the optimization.

    This function uses the dp_model_recomp_heating() function above, but optimization is used
    to determine r0 as the initial starting radius.

    NOTE that the dT_known has to be
    sufficiently smooth for the derivative operation to not be contaminated by noise

    :param fs:                      as the sampling frequency (Hz)
    :param q:                       as the heat input into the soil
    :param k:                       as the thermal conductivity
    :param alpha:
    :param t:                       as the time vector
    :param dT_known:                as the change in temperature
    :param r0_start:                as the initial starting r0
    :param get_gamma4:              True to return gamma4
    :return:
    """
    pass

