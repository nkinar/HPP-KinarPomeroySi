import numpy as np
from scipy.optimize import minimize
from forward_model_nominal import compute_delta_T_dual_heating_cooling, compute_single_probe_forward_model, \
                                compute_single_probe_forward_model_min, compute_delta_T_dual_heating_cooling_q_H, \
                                compute_delta_T_dual_heating_q_H, compute_delta_T_dual_infinite_line_heating
np.seterr(over='ignore')    # sometimes the sum needs to overflow


"""
This file contains the inverse model that is nominally used for the heat pulse probe
"""


def obtain_r_from_curve(q, k, alpha, t, t_heating, known, rstart):
    """
    The input is a fully-processed difference curve.
    The output is a known r found using optimization.
    :param q:
    :param k:
    :param alpha:
    :param t:
     :param t_heating:
    :param known:
    :param rstart:
    :return:            r as the found radius
    """
    def _f_obtain_r_from_curve(x, *args):
        """
        Optimization function to obtain r from the curve.
        NOTE that both the heating and cooling curves are used to obtain the radius.
        :param x:     (r)
        :param args:  (q, k, alpha, t, known)
        :return:
        """
        # variable to find
        r = x[0]
        # known variables
        q = args[0]
        k = args[1]
        alpha = args[2]
        t = args[3]
        t_heating = args[4]
        known = args[5]
        # compute synthetic
        deltaT = compute_delta_T_dual_heating_cooling(q, k, alpha, r, t, t_heating)
        diff = (known - deltaT) ** 2
        dsum = np.sum(diff)
        return dsum
    xstart = (rstart,)
    args = (q, k, alpha, t, t_heating, known)
    res = minimize(_f_obtain_r_from_curve, xstart, args, method='Nelder-Mead')
    rout = res.x[0]
    return rout
# DONE


def obtain_sp_vars_from_curve(q, t, known, kstart):
    """
    Obtain the single probe variables from the curve
    :param q:
    :param kstart:
    :return:
    """
    def _f_obtain_sp(x, *args):
        """
        Curve-fitting for the single-probe
        :param x:
        :param args:
        :return:
        """
        q = args[0]
        t = args[1]
        known = args[2]

        k = x[0]
        B = x[1]
        C = x[2]
        D = x[3]
        deltaT = compute_single_probe_forward_model(q, k, t, B, C, D)
        diff = (known - deltaT) ** 2
        dsum = np.sum(diff)
        return dsum
    B = C = D = 1  # starting values for coefficients
    args = (q, t, known)
    xstart = (kstart, B, C, D)
    res = minimize(_f_obtain_sp, xstart, args, method='Nelder-Mead')
    k_out = res.x[0]
    b_out = res.x[1]
    c_out = res.x[2]
    d_out = res.x[3]
    return k_out, b_out, c_out, d_out
# DONE


def obtain_sp_vars_from_curve_min(q, t, known, kstart, t0start, dstart):
    """
    Obtain the single probe variables from the curve
    :param q:
    :param kstart:
    :param t0start:
    :param dstart:
    :return:
    """
    def _f_obtain_sp_min(x, *args):
        q = args[0]
        t = args[1]
        known = args[2]

        k = x[0]
        t0 = x[1]
        d = x[2]

        deltaT = compute_single_probe_forward_model_min(q, k, t, t0, d)
        diff = (known - deltaT) ** 2
        dsum = np.sum(diff)
        return dsum
    args = (q, t, known)
    xstart = (kstart, t0start, dstart)
    res = minimize(_f_obtain_sp_min, xstart, args, method='Nelder-Mead')
    k_out = res.x[0]
    t0_out = res.x[1]
    d_out = res.x[2]
    return k_out, t0_out, d_out
# DONE


def obtain_k_alpha_from_dual_probe(q, t, t_heating, known, kstart, Hstart, r, full=True):
    """
    Obtain the (k, alpha) from the dual probe model.
    NOTE that if full is True, then the full model (heating and cooling) is used.

    :param q:           as the heat input into the soil
    :param t:           as the time vector
    :param t_heating:   as the time of heating
    :param known:       as the known curve
    :param kstart:      as the starting k
    :param Hstart:      as the starting H
    :param r:           as the known assumed radius
    :param full:        True to use the full curve
    :return:
    """
    k, h = obtain_k_H_from_dual_probe(q, t, t_heating, known, kstart, Hstart, full)
    alpha = r**2 / (4.0*h)
    return k, alpha


def obtain_k_H_from_dual_probe(q, t, t_heating, known, kstart, Hstart, full=True):
    """
    Function to obtain k and H from the heating dual probe model
    :param q:               as the heat input into the soil
    :param t:               as the time vector
    :param known:           as the change in temperature
    :param kstart:          as the starting value for k
    :param Hstart:          as the H value that is used in the equation
    :return:
    """
    if full:
        f = compute_delta_T_dual_heating_cooling_q_H
    else:
        f = compute_delta_T_dual_heating_q_H

    def _f_obtain_q_from_curve(x, *args):
        q = args[0]
        t = args[1]
        t_heat = args[2]
        known = args[3]
        k = x[0]
        H = x[1]
        if full:
            deltaT = f(q, k, H, t, t_heating)
        else:
            deltaT = f(q, k, H, t)
        diff = (known - deltaT)**2
        dsum = np.sum(diff)
        return dsum

    args = (q, t, t_heating, known)
    xstart = (kstart, Hstart)
    res = minimize(_f_obtain_q_from_curve, xstart, args, method='Nelder-Mead')
    k_out = res.x[0]
    h_out = res.x[1]
    return k_out, h_out
# DONE


def inverse_model_dual_probe_variable_radius(q, t, t_heating, r_t, known, kstart, alpha_start, full=True):
    """
    Run the inverse model to obtain a dual probe inverse with variable radius.
    This works for
    1. heating                  full = False
    2. heating and cooling      full = True
    3. cooling only             full = True with only the cooling section of the curve being used

    :param q:                       as the known heat input into the soil
    :param t:                       as the time vector
    :param t_heating:               as the time of heating
    :param r_t:                     as the time-variable r_t
    :param known:                   as the known curve
    :param kstart:                  as the starting k
    :param alpha_start:             as the starting alpha
    :param full:                    True to use the full curve
    :return:
    """
    if full:
        f = compute_delta_T_dual_heating_cooling
    else:
        f = compute_delta_T_dual_infinite_line_heating

    def _f_obtain_alpha_r(x, *args):
        q = args[0]
        t = args[1]
        r = args[2]
        known = args[3]
        k = x[0]
        alpha = x[1]
        if full:
            dT = f(q, k, alpha, r, t, t_heating)
        else:
            dT = f(q, k, alpha, r, t)
        diff = (known - dT)**2
        dsum = np.sum(diff)
        print('dsum = ', dsum)
        return dsum
    args = (q, t, r_t, known)
    xstart = (kstart, alpha_start)
    res = minimize(_f_obtain_alpha_r, xstart, args, method='Nelder-Mead')
    k_out = res.x[0]
    alpha_out = res.x[1]
    return k_out, alpha_out










