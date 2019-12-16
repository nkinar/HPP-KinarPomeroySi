import scipy as sp
import numpy as np
from get_size import length
import sys
from get_alpha import get_alpha
from is_numpy import is_numpy
from ei import ei

"""
This file runs the forward model that is nominally used for the Heat Pulse Probe (HPP)
"""
# CONSTANTS
pi = np.pi


def term1(q, k):
    """
    Compute term 1 for the infinite line source
    :param q:
    :param k:
    :return:
    """
    out = -q / (4*pi*k)
    return out


def term2(r, alpha, t):
    """
    Compute term 2 for the infinite line source
    :param r:
    :param alpha:
    :param t:
    :return:
    """
    out = -(r**2) / (4*alpha*t)
    return out


def term3(r, alpha, t, th):
    """
    Compute term3
    :param r:
    :param alpha:
    :param t:
    :param th:
    :return:
    """
    out = -(r**2) / (4 * alpha * (t - th))
    return out

def compute_delta_T_dual_infinite_line_heating(q, k, alpha, r, t):
    """
    Compute the change in temperature for the dual probe
    :param q:           quantity of heat liberated
    :param k:           thermal conductivity
    :param alpha:       thermal diffusivity
    :param r:           probe radius
    :param t:           time
    :return:
    """
    first = term1(q, k)
    second = term2(r, alpha, t)
    term = ei(second)
    deltaT = first * term
    return deltaT


def compute_delta_T_dual_heating_cooling(q, k, alpha, r, t, t_heating, split=False):
    """
    Compute the synthetic curve using the heating and cooling
    :param q:           quantity of heat liberated
    :param k:           thermal conductivity
    :param alpha:       thermal diffusivity
    :param r:           probe radius
    :param t:           time vector
    :param t_heating:   time at which heating ends
    :return:
    """
    n = length(t)
    out = np.zeros(n)
    first = term1(q, k)
    second = -first
    dt = t[1] - t[0]
    fs = 1.0 / dt
    nn = int(np.ceil(fs*t_heating))
    tvec_heat = t[:nn]
    tvec_cool = t[nn:]
    if is_numpy(r):         # variable r
        rheat = r[:nn]
        rcool = r[nn:]
    else:                   # fixed r
        rheat = r
        rcool = r
    heating = compute_delta_T_dual_infinite_line_heating(q, k, alpha, rheat, tvec_heat)
    ei_term3 = ei(term3(rcool, alpha, tvec_cool, t_heating))
    ei_term2 = ei(term2(rcool, alpha, tvec_cool))
    cooling = second * (ei_term3 - ei_term2)
    out = np.append(heating, cooling)
    if split:
        return out, heating, cooling, first, second, tvec_heat, tvec_cool, rheat, rcool, ei_term3, ei_term2
    return out


def compute_delta_T_dual_heating_cooling_q_H(q, k, H, t, t_heating):
    """
    Compute the synthetic curve using the heating and cooling with only H
    :param q:           quantity of heat liberated
    :param k:           thermal conductivity
    :param H:           H parameter
    :param t:           time vector
    :param t_heating:   time at which heating ends
    :return:
    """
    n = length(t)
    out = np.zeros(n)
    first = term1(q, k)
    for k in range(n):
        tk = t[k]
        if tk > t_heating:
            second = ei(-H/(tk-t_heating)) - ei(-H/tk)
            second = -second
        else:
            second = ei(-H/tk)
        out[k] = first * second
    return out


def compute_delta_T_dual_heating_q_H(q, k, H, t):
    """
    Compute the synthetic curve using the heating and cooling with only H
    :param q:           quantity of heat liberated
    :param k:           thermal conductivity
    :param H:           H parameter
    :param t:           time vector
    :param t_heating:   time at which heating ends
    :return:
    """
    n = length(t)
    out = np.zeros(n)
    first = term1(q, k)
    second = ei(-H / t)
    out = first * second
    return out


def compute_single_probe_forward_model(q, k, t, B, C, D, constrain=True):
    """
    Compute the single probe model for heating
    :param q:
    :param k:
    :param t:
    :param B:
    :param C:
    :param D:
    :return:
    """
    log_t = np.log(t)
    first = -term1(q, k) * log_t
    second = B
    third = (1/t)*(C*log_t + D)
    out = first + second + third
    if constrain:
        out[out < 0] = 0
    return out


def compute_single_probe_forward_model_min(q, k, t, t0, d):
    log_t = np.log(t)
    out = -term1(q, k) * log_t + d
    out[out < 0] = 0
    return out


