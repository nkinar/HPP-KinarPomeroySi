import numpy as np
from numpy import pi
from numerical_derivative import first_derivative_3pt, second_derivative_3pt
from scipy.optimize import minimize
from get_size import length


"""
This code is responsible for computing k using the derivative operation and the signal processing
for the Single Probe (SP).
"""


def sp_model_nominal(q, k, t, b, c, d):
    """
    SP forward model
    :param q:
    :param k:
    :param t:
    :param b:
    :param c:
    :param d:
    :return:
    """
    term0 = (q/(4*pi*k))*np.log(t) + b
    term1 = (1/t)*(c*np.log(t) + d)
    out = term0 + term1
    return out


def sp_model_late_time(q, k, t, b):
    """
    Apply the SP forward model with c = d = 0
    :param q:
    :param k:
    :param t:
    :param b:
    :return:
    """
    term0 = (q/(4*pi*k))*np.log(t) + b
    return term0


def sp_model_nominal_get_b_c_d(dT, t, q, k):
    """
    SP model with known {q, k, b} to obtain {c, d}
    :param dT:      as the change in temperature curve
    :param t:       as the time vector
    :param q:       as the known heat input into the soil
    :param k:       as the thermal conductivity
    :param b:       as the slope of the curve
    :return:        (c, d) as the other parameters of the model
    """
    def _f_obtain_sp_c_d(x, *args):
        q = args[0]
        t = args[1]
        k = args[2]
        known = args[2]
        b = x[0]
        c = x[1]
        d = x[2]
        dT = sp_model_nominal(q, k, t, b, c, d)
        diff = (known - dT) ** 2
        dsum = np.sum(diff)
        return dsum
    args = (q, t, k)
    xstart = (1.0, 1.0, 1.0)
    res = minimize(_f_obtain_sp_c_d, xstart, args, method='Nelder-Mead')
    bout = res.x[0]
    cout = res.x[1]
    dout = res.x[2]
    return bout, cout, dout


def sp_model_late_time_inv_optim(dT, t, q):
    """
    Obtain the SP model for late-time application
    :param dT:
    :param t:
    :param q:
    :return:
    """
    def _f_obtain_sp_late(x, *args):
        q = args[0]
        t = args[1]
        known = args[2]
        k = x[0]
        b = x[1]
        dT_calc = sp_model_late_time(q, k, t, b)
        diff = (known - dT_calc)**2
        dsum = np.sum(diff)
        return dsum
    args = (q, t, dT)
    xstart = (1.0, 1.0)
    res = minimize(_f_obtain_sp_late, xstart, args, method='Nelder-Mead')
    kdet = res.x[0]
    bdet = res.x[1]
    return kdet, bdet


def dT_sp_model_signal(q, k, t, b):
    """
    Be able to run the change in temperature for SP single model
    :param q:
    :param k:
    :param t:
    :param b:
    :return:
    """
    dT = (q / (4 * pi * k)) * (np.log(t) + 2.0) + b
    return dT


def sp_model_signal_inv(dT, t, dt, q, output_model_array=False, cutoff=0.3, filter=True):
    """
    Obtain the thermal conductivity k using the derivative signal processing method.
    NOTE that the trim operation is included as well.

    :param dT:                      as the single probe change curve
    :param t:                       as the time vector
    :param dt:                      as the timestep
    :param q:                       as the heat input to the medium
    :param output_model_array:      True to output the model array
    :param cutoff:                  Set the filter cutoff frequency
    :return:
    k as the thermal conductivity that has been detected by the processing
    b as the b value
    """
    from butterworth_low import butter_lowpass_filter

    def _f_obtain_sp_signal(x, *args):

        q = args[0]
        tt = args[1]
        known = args[2]
        k = x[0]
        b = x[1]

        dT = dT_sp_model_signal(q, k, tt, b)
        diff = (known - dT)**2
        dsum = np.sum(diff)
        return dsum
    dT0 = t * dT                                            # homodyne by t
    tt, dT1 = first_derivative_3pt(dT0, dt, t)              # derivative
    dT1_h = tt * dT1                                        # homodyne by t
    tth, dT2_h = first_derivative_3pt(dT1_h, dt, tt)        # derivative
    if filter:  # filter the signal if required
        fs = 1 / dt
        dT2_h = butter_lowpass_filter(dT2_h, cutoff, fs, order=5, both=True)

    args = (q, tth, dT2_h)
    xstart = (1.0, 1.0)
    res = minimize(_f_obtain_sp_signal, xstart, args, method='Nelder-Mead')
    kdet = res.x[0]
    bdet = res.x[1]
    if not output_model_array:
        return kdet, bdet
    return kdet, bdet, tth, dT2_h


##########################################################################################


def sp_model_late_time_inv(dT, t, dt, q, entire_heating_time=False, return_all=False):
    """
    Obtain the k value from the late-time inverse
    :param dT:                          change in temperature
    :param t:                           time vector
    :param dt:                          time step
    :param q:                           heat input into soil
    :param entire_heating_time:         True to not search for when the curve becomes linear
    :return: (tlinear, kdet)
    tlinear         = time at which the curve becomes linear
    kdet            = determined k
    """
    TADD = 1.0          # time additional to take when the curve is linear (ensures linearity)
    log_t = np.log(t)   # t cannot be zero since we are taking the logarithm
    if not entire_heating_time:
        tt, sd = second_derivative_3pt(dT, dt, log_t)
        nn = length(sd)
        kk = 0
        for k in range(1, nn):      # zero-crossing detector to determine when the curve is linear
            check = sd[k]           # start at the first element
            if check < 0:
                kk = k
                break
        k2 = kk-2
        if k2 < 0:                  # ensure that the element selected is not less than 0
            k2 = 0
        tlinear = t[k2] + TADD      # time at which the curve becomes linear
        fs = 1.0/dt                 # sampling frequency used to compute the cut time
        tlinear_cut_elem = int(np.ceil(tlinear*fs))
        t_tlinear = t[tlinear_cut_elem:]
        dT_tlinear = dT[tlinear_cut_elem:]
    else:  # bypass cutting the curve
        t_tlinear = t
        dT_tlinear = dT
        tlinear = 0

    # take the curve with the late time since this is fit over the late time section of the curve
    t_tlinear0 = t_tlinear
    dT_tlinear0 = dT_tlinear

    kdet, bdet = sp_model_late_time_inv_optim(dT_tlinear0, t_tlinear0, q)
    if not return_all:
        return tlinear, kdet
    return tlinear, kdet, bdet, t_tlinear0, dT_tlinear0


def main():
    pass


if __name__ == '__main__':
    main()

