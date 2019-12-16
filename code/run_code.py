import numpy as np
from forward_model_nominal import compute_delta_T_dual_infinite_line_heating, get_alpha, \
    compute_delta_T_dual_heating_cooling, compute_single_probe_forward_model
from inverse_model_nominal import obtain_r_from_curve, obtain_sp_vars_from_curve
from get_size import length
from numpy import genfromtxt
from scipy.signal import medfilt
# from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.signal import butter, filtfilt
from average_downsample import average_downsample
import sys

# common frontend to place at the top of each
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('ggplot')
mpl.rcParams["mathtext.fontset"] = "stix"

from constants import *


######################################################################################################################


from forward_model_nominal import compute_expi_array
def compute_e1(x):
    return -compute_expi_array(-x)


def plot_ei():
    n = 100
    x = np.linspace(0.01, 2, n)
    y = compute_e1(x)
    left = -1/np.sqrt(x)
    right = 1/np.sqrt(x)

    left0 = -1 / compute_e1(x)
    right0 = 1 / compute_e1(x)

    block = False
    plt.figure()
    plt.plot(x, y)
    plt.plot(x, left)
    plt.plot(x, right)
    plt.show(block)

    block = True
    plt.figure()
    plt.plot(-x, -x)
    plt.plot(-x, left0)
    plt.plot(-x, right0)
    plt.show(block)


def plot_ineq():
    n = 100
    h = 80
    t = np.linspace(5, 10, n)
    th = 5
    tth = t - th
    t = t[1:]
    tth = tth[1:]
    x = np.sqrt(h/tth)
    y = np.sqrt(h/t)
    term = -compute_e1(x) + compute_e1(y)
    left = -1/term
    right = 1/term
    mid = 1/np.sqrt(x) + 1/np.sqrt(y)

    block = True
    plt.figure()
    plt.plot(t, mid)
    plt.plot(t, left)
    plt.plot(t, right)
    plt.show(block)


def check_Eix():
    n = 100
    x = np.linspace(0.01, 2, n)
    left = -1 / np.sqrt(x)
    right = 1 / np.sqrt(x)
    y = -compute_e1(x)

    block = True
    plt.figure()
    plt.plot(x, left)
    plt.plot(x, right)
    plt.plot(x, y)
    plt.show(block)

def plot_ineq_test():
    from ei import ei, e1
    n = 100
    x = np.linspace(1, 0.001, n)
    y = e1(x)
    left = -1.0 / np.sqrt(x)
    right = 1.0 / np.sqrt(x)

    block = True
    plt.figure()
    plt.plot(x, y)
    plt.plot(x, left)
    plt.plot(x, right)
    plt.show(block)


def plot_ineq_check():
    """
    Check the second part of the inequality
    :return:
    """
    from forward_model_nominal import compute_delta_T_dual_heating_cooling
    from get_volumetric_heat_capacity import get_volumetric_heat_capacity
    from gen_time_vec import gen_time_vec
    from get_alpha import get_alpha_k_C
    from numpy import pi
    from ei_inv_cooling import ei_inv_cooling_array

    r = 6e-3
    q = 50.0
    k = 5.0
    theta_o = 0.01
    theta_w = 0.40
    theta_m = 1 - theta_o - theta_w
    c = get_volumetric_heat_capacity(theta_m, theta_o, theta_w)
    alpha = get_alpha_k_C(k, c)
    t_heating = 8
    th = t_heating
    T = 30
    fs = 120
    t = gen_time_vec(fs, T, tstart=0)
    t = t[1:]
    dT_hc = compute_delta_T_dual_heating_cooling(q, k, alpha, r, t, t_heating)
    nn = int(np.ceil(fs*t_heating))
    dT_cool = dT_hc[nn:]
    t_cool = t[nn:]

    term = (4*pi*k) / q
    dT_cool0 = term*dT_cool  # this is the input to find the inverse

    h = (r**2)/(4*alpha)
    x = h / (t_cool-th)
    y = h / t
    z = np.sqrt((t_cool-th)/t_cool)
    high = (t_cool-th)*(((1+z)/z)*(1.0/dT_cool0))
    print('h = ', h)
    aout = ei_inv_cooling_array(dT_cool0, t_cool, th)
    print('done inverse')

    block = False
    plt.figure()
    plt.plot(t, dT_hc)
    plt.plot(t_cool, dT_cool)
    plt.axvline(x=t_heating)
    plt.show(block)

    block = False
    plt.figure()
    plt.plot(t_cool, dT_cool0)
    plt.plot(t_cool, high)
    plt.show(block)

    block = True
    plt.figure()
    plt.plot(t_cool, aout)
    plt.show(block)

################################################################################

from deriv_forward_inv import inverse_first_derivative


def _fp(x, *args):
    # UNKNOWN
    r0 = x[0]               # r0 as the unknown starting radius
    alpha = x[1]            # as the unknown alpha

    # KNOWN
    dT = args[0]            # change in temperature
    k = args[1]             # thermal conductivity
    t = args[2]             # time
    q = args[3]             # q
    t_heating = args[4]     # time of heating
    gamma7 = args[5]        # known gamma7
    dt = args[6]            # timestep

    # compute the dT for comparison
    r = inverse_first_derivative(gamma7, dt, np.log(r0))
    # compute the forward model
    dT_comp = compute_delta_T_dual_heating_cooling(q, k, alpha, r, t, t_heating, split=False)
    # compute the difference
    diff = (dT_comp-dT)**2
    # we want to minimize the sum of squared errors
    dsum = np.sum(diff)
    return dsum




################################################################################


def run_find_hpp_nominal():
    """
    Do experiment to obtain nominal HPP data
    :return:
    """
    from calibration_processing import initial_compute_hpp
    from comparisons import compute_rmse, compute_mb
    from inverse_model_nominal import obtain_sp_vars_from_curve_min, compute_single_probe_forward_model_min
    # path = '../hpp-data/hpp-formatted-data/calibration/July-12-2017/'
    path = '../hpp-data/hpp-formatted-data/sand/July-20-2017/'
    # path = '../hpp-data/hpp-formatted-data/sand/July-21-2017/'
    # path = '../hpp-data/hpp-formatted-data/sand/July-26-2017/'
    # path = '../hpp-data/hpp-formatted-data/peat/July-26-2017/'
    # path = '../hpp-data/hpp-formatted-data/peat/July-27-2017/'
    # path = '../hpp-data/hpp-formatted-data/peat/July-30-2017/'
    # path = '../hpp-data/hpp-formatted-data/peat/July-30-2017/'
    # path = '../hpp-data/hpp-formatted-data/peat/July-31-2017/'
    start_string = 'sand-'
    # start_string = 'cal-'
    # start_string = 'peat-'
    downsample = 'dp'
    filt = 'both'
    filt_cutoff = 10            # cutoff filter frequency (Hz)
    fs = 120                    # sample rate of device (Hz)
    fs_down = 12                # downsample frequency (Hz)
    t_cold = 1                  # time of cold sampling at beginning (s)
    additional_text = '-rep1'
    r0 = 6e-3                   # starting radius (m)

    # Tikhonov regularization parameter
    epsilon = 0.0

    #############################################################################

    q = 45
    t_heat = 8
    t_total = 3 * SECONDS_IN_MIN

    # q = 10
    # t_heat = 10 * SECONDS_IN_MIN
    # t_total = 12 * SECONDS_IN_MIN

    ################################################

    # LOAD IN THE COEFFICIENTS USED FOR CALIBRATION AND OBTAIN r_nominal
    # NOTE that the calibrated r0 is stored in mm, so we need to convert to m to be able to use the data
    from obtain_load_r0 import SaveLoadR0
    sr = SaveLoadR0()
    sr.load_r0(CONST_PATH + CAL_R0_FILENAME + SAVEZ_EXT)
    r_nominal = sr.get_r0(q) * 1.0e-3
    # r_nominal = 6e-3

    # We always use assumed q and use the step detector q
    use_assumed_q = False
    use_step_detector_q = True

    I, Rw, T1, T2, Vout, dV, dac_code, delta_T1, delta_T1_trim, delta_T1_trim_heating, \
    delta_T2, delta_T2_trim, num, qav, qprime, step_q_tfirst, step_q_tsecond, t1, t1_trim, \
    t1_trim_heating, t2, t2_trim, ypk, t, rmse_q_calc, mb_q_calc, pd_q_calc, t2_trim_heating, \
    delta_T2_trim_heating = \
        initial_compute_hpp(start_string, additional_text,
        downsample, filt, filt_cutoff, fs, fs_down, path, q, t_cold, t_heat, t_total, use_assumed_q,
        use_step_detector_q)

    dt_sp = t1_trim_heating[1] - t1_trim_heating[0]
    dt_dp = t2_trim_heating[1] - t2_trim_heating[0]

    from sp_model_signal import sp_model_signal_inv
    from curve_time_trim import curve_time_trim
    from sp_model_signal import sp_model_late_time_inv
    from inverse_model_nominal import obtain_k_alpha_from_dual_probe
    from inverse_model_dp_radius import InverseWarmupCurveDP
    from get_theta_rho import get_theta_rho

    # heating only for SP
    t1_trim_heating1, delta_T1_trim_heating1 = curve_time_trim(t1_trim_heating, delta_T1_trim_heating)
    # heating only for DP
    t2_trim_heating1, delta_T2_trim_heating1 = curve_time_trim(t2_trim_heating, delta_T2_trim_heating)
    # heating and cooling for DP
    t2_trim1, delta_T2_trim1 = curve_time_trim(t2_trim, delta_T2_trim)

    # SP SIGNAL PROCESSING INVERSE MODEL
    kdet_sig, bdet = sp_model_signal_inv(delta_T1_trim_heating1, t1_trim_heating1, dt_sp, q, output_model_array=False)
    # SP inverse model to obtain k via a late-time model
    tlinear, kout_sp_late_time = sp_model_late_time_inv(delta_T1_trim_heating1, t1_trim_heating1, dt_sp, qav)

    # Calculate the time at which the heat pulse goes linear relative to the start of the experiment.
    # This also includes the cold time at the beginning of the experiment (used for graphing).
    tlinear_from_begin = t_cold + tlinear

    # SP FORWARD MODELS
    kd_signal, bd_signal, cd_signal, dd_signal = obtain_sp_vars_from_curve(qav,
                                                                           t1_trim_heating1,
                                                                           delta_T1_trim_heating1,
                                                                           kout_sp_late_time)
    kd_lt, bd_lt, cd_lt, dd_lt = obtain_sp_vars_from_curve(qav, t1_trim_heating1, delta_T1_trim_heating1,
                                                           kout_sp_late_time)

    # SP SYNTHETIC MODELS
    delta_T_single_synth_signal = compute_single_probe_forward_model(qav, kd_signal, t1_trim_heating1, bd_signal,
                                                              cd_signal, dd_signal)
    delta_T_single_synth_late_time = compute_single_probe_forward_model(qav, kd_lt, t1_trim_heating1, bd_lt, cd_lt,
                                                                        dd_lt)

    # OBTAIN {k, alpha} FROM THE NOMINAL DP MODEL (HEATING AND COOLING)
    # NOTE THAT THIS REQUIRES AN ESTIMATE OF THE PROBE SPACING RADIUS as r_nominal READ FROM A COEFFICIENT FILE
    kstart = kd_signal
    Hstart = 80
    k_nom_heat, alpha_nom_heat = obtain_k_alpha_from_dual_probe(q, t2_trim_heating1, t_heat, delta_T2_trim_heating1,
                                                                kstart, Hstart, r_nominal, full=False)  # heat only
    k_nom_heat_and_cool, alpha_nom_heat_and_cool = obtain_k_alpha_from_dual_probe(qav, t2_trim1, t_heat, delta_T2_trim1,
                                                                kstart, Hstart, r_nominal, full=True)  # heat and cool
    # OBTAIN {r(t), alpha} FROM THE SIGNAL PROCESSING MODEL
    # r(t) = r_t_sig as the time-variable radius
    InvDP = InverseWarmupCurveDP()
    typ = 'iterative'
    from scipy.interpolate import spline
    from scipy.interpolate import UnivariateSpline
    # s = spline(t2_trim_heating1, delta_T2_trim_heating1, t2_trim_heating1, order=5, kind='smoothest')
    spl = UnivariateSpline(t2_trim_heating1, delta_T2_trim_heating1)
    s = spl(t2_trim_heating1)   # spline (apply to inverse warmup curve)
    from deriv_forward_inv import forward_first_derivative
    sderiv = forward_first_derivative(s, 1/fs_down)

    # run a zero crossing detector to determine when sderiv is zero and then cut at that point
    nn = length(sderiv)
    ncut = 0
    for k in range(nn):
        if sderiv[k] > 0:
            ncut = k
            break
    ncut += 1
    print('ncut = ', ncut)
    tt_ncut = t2_trim_heating1[ncut+1:]
    sderiv_ncut = sderiv[ncut:]
    tt_cut = tt_ncut[0]

    # block = True
    # plt.figure()
    # # plt.plot(tt_ncut, sderiv_ncut)
    # plt.plot(t2_trim_heating1[1:], sderiv)
    # plt.axvline(x=tt_cut)
    # plt.show(block)

    scut = s[ncut:]
    tcut = t2_trim_heating1[ncut:]

    # scut = s
    # tcut = t2_trim_heating1
    # compute from the inverse warmup curve
    r0 = 6e-3
    alpha_sig, r_t_sig = InvDP.inverse_warmup_curve_dp(kdet_sig, scut, qav, tcut,
                                                       r_nominal,
                                                       typ, lowpass=False, epsilon=0.0)

    # average over r_t_sig
    from integrate_average import integrate_average
    r_av = integrate_average(r_t_sig, tcut, 'simpson')
    # r_av = np.mean(r_t_sig)
    print('r_av = ', r_av)

    #############################################################################
    # OBTAIN {k, alpha} FROM THE NOMINAL DP MODEL (HEATING AND COOLING)
    # USE THE r_av FROM SIGNAL PROCESSING
    kstart = kd_signal
    Hstart = 80
    # r_av = r_nominal
    print('r_nominal = ', r_nominal)
    k_sig_heat, alpha_sig_heat = obtain_k_alpha_from_dual_probe(q, t2_trim_heating1, t_heat, delta_T2_trim_heating1,
                                                                kstart, Hstart, r_av, full=False)   # heat only
    k_sig_heat_and_cool, alpha_sig_heat_and_cool = obtain_k_alpha_from_dual_probe(qav, t2_trim1, t_heat, delta_T2_trim1,
                                                                kstart, Hstart, r_av, full=True)    # heat and cool
    #############################################################################


    # DP SYNTHETICS NOMINAL
    deltaT_dual_synth_nominal_heat_cool = compute_delta_T_dual_heating_cooling(qav, k_nom_heat, alpha_nom_heat,
                                                                               r_nominal, t2_trim1, t_heat)
    deltaT_dual_synth_nominal_heating = compute_delta_T_dual_infinite_line_heating(qav, k_nom_heat_and_cool,
                                                                                   alpha_nom_heat_and_cool,
                                                                                   r_nominal, t2_trim_heating1)

    # DP SYNTHETIC WITH RADIUS DETERMINED FROM SIGNAL PROCESSING
    deltaT_dual_syth_variable_radius_heating = compute_delta_T_dual_infinite_line_heating(qav, kdet_sig,
                                                                                         alpha_sig,
                                                                                         r_av,
                                                                                         t2_trim_heating1)
    # set the organic content
    theta_o = 9.2e-3

    # COMPUTATIONS
    # {theta_w, rho} for DP (heating)
    theta_w_heat, rho_heat = get_theta_rho(k_nom_heat, alpha_nom_heat, theta_o)
    # {theta_w, rho} for DP (heating and cooling)
    theta_w_heat_cool, rho_heat_cool = get_theta_rho(k_nom_heat_and_cool, alpha_nom_heat_and_cool, theta_o)

    # SIGNAL PROCESSING (HEAT AND COOL)
    theta_w_heat_sig, rho_heat_sig = get_theta_rho(kdet_sig, alpha_sig, theta_o)
    theta_w_heat_cool_sig, rho_heat_cool_sig = get_theta_rho(kdet_sig, alpha_sig, theta_o)

    # DUAL PROBE
    # RMSE, MB asnd PD for DP model actual vs synthetic (heating)
    # RMSE, MB and PD for DP model actual vs synthetic (heating and cooling)
    # RMSE, MB and PD for DP model actual vs synthetic variable radius (heating only) SIGNAL PROCESSING
    # SINGLE PROBE
    # RMSE, MB and PD for SP model actual vs synthetic (late-time)
    # RMSE, MB and PD for SP model actual vs synthetic (signal processing)

    # from butterworth_low import butter_lowpass_filter
    # y = butter_lowpass_filter(delta_T2_trim_heating1, 1, fs, order=5, both=True)

    # block = False
    # plt.figure()
    # plt.plot(t2_trim_heating1, y)
    # plt.xlabel('Time (s)')
    # plt.ylabel('delta T2 heating')
    # plt.legend()
    # plt.show(block)

    print('OUTPUTS:')
    print('kdet_sig (signal processing) = ', kdet_sig)
    print('tlinear (late-time) = ', tlinear)
    print('kout_sp_late_time (late-time) = ', kout_sp_late_time)
    print('')
    print('k_nom_heat = ', k_nom_heat)
    print('alpha_nom_heat = ', alpha_nom_heat)
    print('')
    print('k_nom_heat_and_cool = ', k_nom_heat_and_cool)
    print('alpha_nom_heat_and_cool = ', alpha_nom_heat_and_cool)
    print('')
    print('alpha_sig (from signal processing with variable radius) = ', alpha_sig)
    print('THETA AND RHO')
    print('-----')
    print('theta_w_heat = ', theta_w_heat)
    print('rho_heat = ', rho_heat)
    print('-----')
    print('theta_w_heat_cool = ', theta_w_heat_cool)
    print('rho_heat_cool = ', rho_heat_cool)

    print('----SIGNAL PROCESSING----')
    print('theta_w_heat_sig = ', theta_w_heat_sig)
    print('rho_heat_sig = ', rho_heat_sig)
    print('-----')
    print('---HEAT AND COOL:---')
    print('theta_w_heat_cool_sig = ', theta_w_heat_cool_sig)
    print('rho_heat_cool_sig = ', rho_heat_cool_sig)
    print('--------------------------')

    block = False
    plt.figure()
    plt.plot(t1_trim_heating1, delta_T1_trim_heating1, label='data')
    plt.plot(t1_trim_heating1, delta_T_single_synth_signal, label='signal')
    plt.plot(t1_trim_heating1, delta_T_single_synth_late_time, label='late time')
    plt.xlabel('Time (s)')
    plt.ylabel('delta T1 heating SP')
    plt.legend()
    plt.show(block)

    # plt.plot(t2_trim1, delta_T2_trim1, label='data heating and cooling')
    # plt.plot(t2_trim_heating1, deltaT_dual_synth_nominal_heating, label='heating')
    # plt.plot(t2_trim1, deltaT_dual_synth_nominal_heat_cool, label='heating and cooling')
    block = False
    plt.figure()
    plt.plot(t2_trim_heating1, delta_T2_trim_heating1, label='data heating')
    # plt.plot(t2_trim_heating1, deltaT_dual_syth_variable_radius_heating, label='heating (variable radius)')
    plt.plot(t2_trim_heating1, s, label='heating (variable radius) smooth')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('delta T2 heating DP')
    plt.show(block)

    block = True
    plt.figure()
    plt.plot(tcut, r_t_sig, label='r(t)')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('r(t)')
    plt.show(block)


def plot_ei_test():
    from ei import ei
    from get_alpha import get_alpha_k_C
    from get_volumetric_heat_capacity import get_volumetric_heat_capacity
    from gen_time_vec import gen_time_vec

    q = 50
    k = 4
    r0 = 6.7e-3
    r1 = 15e-3
    theta_o = 1.0e-3
    theta_w = 0.49
    theta_m = 1 - theta_o - theta_w
    C = get_volumetric_heat_capacity(theta_m, theta_o, theta_w)
    fs = 12
    dt = 1/fs
    T = 30
    t = gen_time_vec(fs, T)
    t = t[1:]   # cannot start with zero
    n = length(t)
    r = np.linspace(r0, r1, n)
    alpha = get_alpha_k_C(k, C)
    print('alpha = ', alpha)
    t_heating = 8
    th = t_heating
    dT, heating, cooling, first, second, tvec_heat, tvec_cool, rheat, rcool, ei_term3, ei_term2 = \
        compute_delta_T_dual_heating_cooling(q, k, alpha, r, t, t_heating, split=True)

    # check the recomposition
    from dp_model_recomp import dp_model_recomp_heating_cooling, dp_model_recomp_heating
    r_t = dp_model_recomp_heating_cooling(fs, q, k, t, th, dT, r0)
    r_t_heat = dp_model_recomp_heating(fs, q, k, tvec_heat, heating, r0)

    block = True
    plt.figure()
    plt.plot(t, r_t)
    plt.plot(t, r)
    plt.plot(tvec_heat, r_t_heat)
    plt.show(block)

    # dT_known = dT
    # from numpy import pi
    # from dp_model_recomp import ei_inv_neg_array, ei_inv_cooling
    # from deriv_forward_inv import forward_first_derivative, inverse_first_derivative
    # out = np.zeros(n)
    # term = -(4*pi*k)/q
    # dt = 1 / fs
    # nn = int(np.ceil(fs*th))
    # tvec_heat = t[:nn]
    # tvec_cool = t[nn:]
    # dT_known_heat = dT_known[:nn]
    # dT_known_cool = dT_known[nn:]
    # # HEATING
    # gamma2_heat = term*dT_known_heat
    # gamma3_heat = ei_inv_neg_array(gamma2_heat) * tvec_heat
    # gamma4_heat = -gamma3_heat
    # # COOLING
    # dT_cool_strip = -term*dT_known_cool
    # gamma4_cool = ei_inv_cooling(dT_cool_strip, th, tvec_cool)
    # gamma4 = np.concatenate((gamma4_heat, gamma4_cool))   # AFTER HERE
    # gamma5 = np.sqrt(gamma4)
    # gamma6 = np.log(gamma5)
    # gamma7 = forward_first_derivative(gamma6, dt)
    #
    # gamma7_comp = forward_first_derivative(np.log(r), dt)
    # gamma7_diff = gamma7 - gamma7_comp
    #
    # # check the inverse
    # r0_wanted = 6e-3
    # log_r_inv = inverse_first_derivative(gamma7, dt, np.log(r0_wanted))
    # r_inv = np.exp(log_r_inv)
    #
    #
    # block = False
    # plt.figure()
    # plt.plot(t[1:], gamma7)
    # plt.plot(t[1:], gamma7_comp)
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(t[1:], gamma7_diff)
    # plt.show(block)

    # block = True
    # plt.figure()
    # plt.plot(t[1:], linear_test_vec_log_deriv)
    # plt.plot(t[1:], gamma7)
    # plt.show(block)

    # block = True
    # plt.figure()
    # plt.plot(t, r_inv)
    # plt.show(block)


    from dp_model_recomp import dp_model_recomp_heating_cooling
    # dp_model_recomp_heating_cooling(fs, q, k, t, th, dT, find_r0=True, r0=None)
    # diff = r - r_t
    #
    # block = False
    # plt.figure()
    # plt.plot(t, r_t)
    # plt.plot(t, r)
    # plt.show(block)
    #
    # block = True
    # plt.figure()
    # plt.plot(t, diff)
    # plt.show(block)
    #
    # _fp

    # from deriv_forward_inv import forward_first_derivative
    # gamma6 = np.log(r)
    # gamma7 = forward_first_derivative(gamma6, dt)
    #
    # plt.figure()

    # cool0 = cooling / second
    # from dp_model_recomp import dp_model_recomp_heating_cooling
    # r_t_recomp = dp_model_recomp_heating_cooling(fs, q, k, t, th, dT, r0)
    # diff = r_t_recomp - r
    #
    # block = False
    # plt.figure()
    # plt.plot(t, r_t_recomp)
    # plt.plot(t, r)
    # plt.show(block)
    #
    # block = True
    # plt.figure()
    # plt.plot(t, diff)
    # plt.show(block)

    # tc = tvec_cool
    # tth = tc - th
    # h = rcool**2/(4*alpha)
    # x = h/tth
    # y = h/tc
    # fxy = ei(-x) - ei(-y)
    #
    # from scipy.optimize import minimize
    # np.seterr(invalid='ignore')
    #
    # def _func(x, *args):
    #     hh = x[0]
    #     known = args[0]
    #     th_known = args[1]
    #     t_known = args[2]
    #     xx = hh / (t_known-th_known)
    #     yy = hh / t_known
    #     first = ei(-xx)
    #     second = ei(-yy)
    #     calc = first - second   # this operation can cause overflow or underflow
    #     out = known-calc
    #     if out < 0:
    #         out = -out
    #     return out
    #
    # nn = length(fxy)
    # recon = np.zeros(nn)
    # hstart = 80
    # xstart = (hstart,)
    # for k in range(nn):
    #     args = (fxy[k], th, tc[k])
    #     res = minimize(_func, xstart, args, method='Nelder-Mead', tol=1.0e-30)  # NOTE THE TOLERANCE
    #     recon[k] = res.x[0]
    # diff = recon-h
    #
    # gamma5 = recon
    # gamma6 = np.log(np.sqrt(gamma5))
    # gamma6_comp = np.log(rcool) - np.log(2*np.sqrt(alpha))
    # gamma6_diff = gamma6 - gamma6_comp
    #
    # from deriv_forward_inv import forward_first_derivative, inverse_first_derivative
    # dt = t[1] - t[0]
    # gamma7 = forward_first_derivative(gamma6, dt)
    # gamma7_comp = forward_first_derivative(np.log(rcool), dt)
    # gamma7_diff = gamma7 - gamma7_comp
    #
    # gamma8 = inverse_first_derivative(gamma7, dt, np.log(rcool[0]))
    # gamma8_comp = np.log(rcool)
    # gamma8_diff = gamma8 - gamma8_comp
    #
    # r_t_recomp = np.exp(gamma8)
    # r_t_diff = r_t_recomp - rcool
    #
    # #########################################################################
    #
    # block = False
    # plt.figure()
    # # plt.plot(tvec_cool, cool0)
    # plt.plot(tvec_cool, fxy)
    # plt.plot(tvec_cool, ei_term3)
    # plt.plot(tvec_cool, ei_term2)
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(tvec_cool, recon)
    # plt.plot(tvec_cool, h)
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(tvec_cool, diff)
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(tvec_cool, gamma6)
    # plt.plot(tvec_cool, gamma6_comp)
    # plt.title('gamma6')
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(gamma6_diff)
    # plt.title('gamma6_diff')
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(tvec_cool[1:], gamma7)
    # plt.plot(tvec_cool[1:], gamma7_comp)
    # plt.title('gamma7')
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(gamma7_diff)
    # plt.title('gamma7_diff')
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(tvec_cool, gamma8)
    # plt.plot(tvec_cool, gamma8_comp)
    # plt.title('gamma8')
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(gamma8_diff)
    # plt.title('gamma8_diff')
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(tvec_cool, r_t_recomp)
    # plt.plot(tvec_cool, rcool)
    # plt.title('r(t)')
    # plt.show(block)
    #
    # block = True
    # plt.figure()
    # plt.plot(r_t_diff)
    # plt.title('r(t) diff')
    # plt.show(block)

    # import mpmath as mp
    # from mpmath import mpf
    # mp.eps = 100000
    # from gen_time_vec import gen_time_vec
    # # h = 3750000
    # h2 = 80
    # h1 = 4
    # nn = 100
    # h = np.linspace(h1, h2, nn)
    # th = 8
    # fs = 120
    # dt = 1/fs
    # T = 60
    # # t = gen_time_vec(th+dt, fs, T)
    # t = 9
    # tth = t - th
    # x = h / tth
    # y = h / t
    # known = ei(-x) - ei(-y)
    # first = ei(-x)
    # second = ei(-y)
    # # known = (np.exp(h / t) / t) - (np.exp(h / tth) / tth)
    #
    # block = True
    # plt.figure()
    # plt.plot(h, known)
    # plt.show(block)

    # def _func(x, *args):
    #     hh = x[0]
    #     k = args[0]
    #     xx = hh / tth
    #     yy = hh / t
    #     calc = mp.ei(-xx) - mp.ei(-yy)
    #     print(calc)
    #     # calc = (np.exp(hh/t)/t) - (np.exp(hh/tth)/tth)
    #     out = mpf(k)-mpf(calc)
    #     if out < 0:
    #         out = -out
    #     return out
    #
    # hstart = 80
    # args = (known,)
    # xstart = (hstart,)
    # from scipy.optimize import minimize
    # res = minimize(_func, xstart, args, method='Nelder-Mead')
    # print(res.x[0])

    # block = True
    # plt.figure()
    # plt.plot(t, known)
    # plt.show(block)


###################################################################################################


def test_function_magic_bound():
    from MagicBounds import MagicBounds, PARAMETER_NAME, MIN_VAL, MAX_VAL

    def func(r, k, c, rho):
        out = (r**2 * rho * c) / (4.0*k)
        return out

    test_param = [
        {
            PARAMETER_NAME: 'r',
            MIN_VAL: 1.0e-3,
            MAX_VAL: 2e-2
        },
        {
            PARAMETER_NAME: 'k',
            MIN_VAL: 0.4,
            MAX_VAL: 20
        },
        {
            PARAMETER_NAME: 'c',
            MIN_VAL: 1.0e6,
            MAX_VAL: 5.0e6
        },
        {
            PARAMETER_NAME: 'rho',
            MIN_VAL: 900,
            MAX_VAL: 3000
        },
    ]
    mb = MagicBounds(func, test_param)
    h_min_val, h_max_val = mb.run_iterate(quiet=False)
    from ei import ei
    from get_eps import get_eps
    eps = get_eps()

    def test_func(h, t, th):
        tth = t - th
        tth[tth == 0] = eps
        xx = h / tth
        yy = h / t
        first = ei(-xx)
        second = ei(-yy)
        calc = first - second

    from gen_time_vec import gen_time_vec
    fs = 120
    T = 500
    th = 8
    tstart = th + (1/fs)
    t = gen_time_vec(fs, T, tstart=tstart)
    term = (t-th)/t

    block = True
    plt.figure()
    plt.plot(t, term)
    plt.show(block)
# DONE


#############################################################################################################


def run_find_hpp_nominal_test_sig():
    """
    Do experiment to obtain nominal HPP data
    :return:
    """
    from calibration_processing import initial_compute_hpp
    from comparisons import compute_rmse, compute_mb
    from inverse_model_nominal import obtain_sp_vars_from_curve_min, compute_single_probe_forward_model_min
    # path = '../hpp-data/hpp-formatted-data/calibration/July-12-2017/'
    path = '../hpp-data/hpp-formatted-data/sand/July-20-2017/'
    # path = '../hpp-data/hpp-formatted-data/sand/July-21-2017/'
    # path = '../hpp-data/hpp-formatted-data/sand/July-26-2017/'
    # path = '../hpp-data/hpp-formatted-data/peat/July-26-2017/'
    # path = '../hpp-data/hpp-formatted-data/peat/July-27-2017/'
    # path = '../hpp-data/hpp-formatted-data/peat/July-30-2017/'
    # path = '../hpp-data/hpp-formatted-data/peat/July-31-2017/'
    start_string = 'sand-'
    # start_string = 'cal-'
    # start_string = 'sand-'
    # start_string = 'peat-'
    downsample = 'dp'
    filt = 'both'
    filt_cutoff = 10                # cutoff filter frequency (Hz)
    fs = 120                        # sample rate of device (Hz)
    fs_down = 12                    # downsample frequency (Hz)
    t_cold = 1                      # time of cold sampling at beginning (s)
    additional_text = '-rep1'
    r0 = 6e-3                       # starting radius

    # sand
    theta_o = 9.2e-3                # organic content
    theta_m = 0.59                  # mineral content

    # Cm_set = 1.1e6
    # Co_set = 1.0e6
    Cm_set = None
    Co_set = None

    # peat
    # theta_o = 0.49
    # theta_m = 0.01

    #############################################################################

    q = 45
    t_heat = 8
    t_total = 3 * SECONDS_IN_MIN

    # q = 45
    # t_heat = 11
    # t_total = 3 * SECONDS_IN_MIN

    # q = 55
    # t_heat = 20
    # t_total = 3 * SECONDS_IN_MIN

    # q = 55
    # t_heat = 89
    # t_total = 3 * SECONDS_IN_MIN

    # q = 10
    # t_heat = 10 * SECONDS_IN_MIN
    # t_total = 12 * SECONDS_IN_MIN

    # q = 20
    # t_heat = 89
    # t_total = 3 * SECONDS_IN_MIN

    ################################################

    # LOAD IN THE COEFFICIENTS USED FOR CALIBRATION AND OBTAIN r_nominal
    # NOTE that the calibrated r0 is stored in mm, so we need to convert to m to be able to use the data
    from obtain_load_r0 import SaveLoadR0
    sr = SaveLoadR0()
    sr.load_r0(CONST_PATH + CAL_R0_FILENAME + SAVEZ_EXT)
    r_nominal = sr.get_r0(q) * 1.0e-3

    # We always use the assumed q and the step detector q
    use_assumed_q = False
    use_step_detector_q = True

    I, Rw, T1, T2, Vout, dV, dac_code, delta_T1, delta_T1_trim, delta_T1_trim_heating, \
    delta_T2, delta_T2_trim, num, qav, qprime, step_q_tfirst, step_q_tsecond, t1, t1_trim, \
    t1_trim_heating, t2, t2_trim, ypk, t, rmse_q_calc, mb_q_calc, pd_q_calc, t2_trim_heating, \
    delta_T2_trim_heating = \
        initial_compute_hpp(start_string, additional_text,
        downsample, filt, filt_cutoff, fs, fs_down, path, q, t_cold, t_heat, t_total, use_assumed_q,
        use_step_detector_q)

    dt_sp = t1_trim_heating[1] - t1_trim_heating[0]
    dt_dp = t2_trim_heating[1] - t2_trim_heating[0]

    from sp_model_signal import sp_model_signal_inv
    from curve_time_trim import curve_time_trim
    from sp_model_signal import sp_model_late_time_inv
    from inverse_model_nominal import obtain_k_alpha_from_dual_probe
    from inverse_model_dp_radius import InverseWarmupCurveDP
    from get_theta_rho import get_theta_rho

    # heating only for SP
    t1_trim_heating1, delta_T1_trim_heating1 = curve_time_trim(t1_trim_heating, delta_T1_trim_heating)
    # heating only for DP
    t2_trim_heating1, delta_T2_trim_heating1 = curve_time_trim(t2_trim_heating, delta_T2_trim_heating)
    # heating and cooling for DP
    t2_trim1, delta_T2_trim1 = curve_time_trim(t2_trim, delta_T2_trim)

    # SP SIGNAL PROCESSING INVERSE MODEL
    kdet_sig, bdet = sp_model_signal_inv(delta_T1_trim_heating1, t1_trim_heating1, dt_sp, q, output_model_array=False)

    th = t_heat
    kstart = 5
    Hstart = 80
    k_det, alpha_det = obtain_k_alpha_from_dual_probe(qav, t2_trim1, th, delta_T2_trim1, kstart, Hstart,
                                                      r_nominal, full=True)
    dT_synth = compute_delta_T_dual_heating_cooling(qav, k_det, alpha_det, r_nominal, t2_trim1, th, split=False)

    # from lmfit.models import SkewedGaussianModel
    # model = SkewedGaussianModel()
    # params = model.make_params(amplitude=10, center=0, sigma=1, gamma=0)
    # result = model.fit(delta_T2_trim1, params, x=t2_trim1)
    # y = result.best_fit

    from get_max import get_max, get_min
    idx, m = get_max(dT_synth)
    print('dt_dp = ', dt_dp)
    idx_time = dt_dp * idx
    print('idx = ', idx)
    print('idx_time = ', idx_time)

    from deriv_forward_inv import forward_first_derivative
    dT_deriv = forward_first_derivative(dT_synth, dt_dp)
    from numerical_derivative import second_derivative_3pt
    tt, dT_deriv2 = second_derivative_3pt(dT_synth, dt_dp, t2_trim1)

    # from find_plateau_curvature import obtain_curvature_1d
    # cc = obtain_curvature_1d(dT_synth, 5, 3, dt_dp)
    idx_peak, md_peak = get_max(dT_synth)
    idxd, md = get_min(dT_deriv2)
    # if idxd > idx_peak:
    #     idxd_found = idxd

    idxd_found = idx_peak
    idxd_time = dt_dp * idxd_found
    # print('idxd_time = ', idxd_time)
    # print('idxd_time_heat = ', idxd_time + th)

    # set the trim time at the end
    # end_trim_time = 40      # only for 21 July
    end_trim_time = 0
    end_cut_num = int(np.floor(fs_down * end_trim_time))

    # Cut the sequence at the beginning
    # cut_time = idxd_time
    # cut_time = idxd_time
    cut_time = idxd_time
    print('cut_time = ', cut_time)
    cut_num = int(np.floor(fs_down * cut_time))
    if end_cut_num > 0:
        t2_trim1_cut = t2_trim1[cut_num:-end_cut_num]
        delta_T2_trim1_cut = delta_T2_trim1[cut_num:-end_cut_num]
    else:
        t2_trim1_cut = t2_trim1[cut_num:]
        delta_T2_trim1_cut = delta_T2_trim1[cut_num:]

    block = False
    plt.figure()
    plt.plot(t2_trim1, delta_T2_trim1)
    plt.plot(t2_trim1, dT_synth)
    # plt.plot(t2_trim1, y)
    plt.plot(t2_trim1[1:], dT_deriv)
    plt.plot(tt, dT_deriv2)
    # plt.plot(t2_trim1, cc)
    plt.axvline(cut_time)
    plt.title('DP heating and cooling')
    plt.show(block)

    # obtain the radius
    from dp_model_recomp import dp_model_recomp_heating, dp_model_recomp_heating_cooling, \
        dp_model_recomp_heating_cooling_trim
    # r_t_heating = dp_model_recomp_heating(fs_down, q, kdet_sig, t2_trim_heating1_cut, delta_T2_trim_heating1_cut,
    #                                       r_nominal)
    r_t_heating_cooling, gamma4 = dp_model_recomp_heating_cooling_trim(fs_down, qav, kdet_sig, t2_trim1_cut, th,
                                                                    delta_T2_trim1_cut, r_nominal, get_gamma4=True)

    # tcut_heat = 7.0
    # cut_num_heat = int(np.floor(fs_down * tcut_heat))
    # r_t_heating, gamma4_heating = dp_model_recomp_heating(fs_down, qav, kdet_sig, t2_trim_heating1[1:],
    #                                                       delta_T2_trim_heating1[1:], r_nominal, get_gamma4=True)

    # compute alpha
    alpha = (r_t_heating_cooling**2) / (4.0*gamma4)
    alpha_sig = np.mean(alpha)

    print('k_det = ', k_det)
    print('alpha_det = ', alpha_det)

    print('kdet_sig = ', kdet_sig)
    print('alpha_sig = ', alpha_sig)

    print('alpha = ', alpha_sig)

    theta_w_var, rho_nom_var = get_theta_rho(k_det, alpha_det, theta_o, theta_m, Cm_set, Co_set)
    theta_w_sig, rho_w_sig = get_theta_rho(kdet_sig, alpha_sig, theta_o, theta_m, Cm_set, Co_set)

    print('------------------------------')
    print('theta_w_var = ', theta_w_var)
    print('rho_nom_var = ', rho_nom_var)
    print('------------------------------')
    print('theta_w_sig = ', theta_w_sig)
    print('rho_w_sig = ', rho_w_sig)
    print('------------------------------')

    # TEST THE INVERSE
    print('testing the inverse')
    from recomp_theta_rho import recomp_theta_rho_additional
    theta_w_var, rho_nom_var, theta_w_sig, rho_w_sig, \
    t2_trim1, delta_T2_trim1, dT_synth_nom, dT_synth_sig, \
    cut_time, t2_trim1_cut, r_t_heating_cooling, \
    alpha_det, alpha_sig, k_det, kdet_sig, \
    t1_trim_heating1, delta_T1_trim_heating1, dT_synth_sp, \
    delta_T1_trim, t1_trim, qav = \
    recomp_theta_rho_additional(path, start_string, additional_text, q, t_heat, t_total, end_trim_time)
    print('done testing')

    print('------------------------------')
    print('theta_w_var = ', theta_w_var)
    print('rho_nom_var = ', rho_nom_var)
    print('------------------------------')
    print('theta_w_sig = ', theta_w_sig)
    print('rho_w_sig = ', rho_w_sig)
    print('------------------------------')

#######################################################################################################

from dp_model_recomp import dp_model_recomp_shared
from deriv_forward_inv import forward_first_derivative
from ei import ei
from get_alpha import get_alpha_k_C
from get_volumetric_heat_capacity import get_volumetric_heat_capacity
from gen_time_vec import gen_time_vec
from dp_model_recomp import dp_model_recomp_heating
from deriv_forward_inv import forward_first_derivative
from ei import ei_inv_neg_array
pi = np.pi
from scipy.optimize import minimize


def _recomp(x, *args):
    r0 = x[0]           # UNKNOWN
    k = args[0]         # known
    t = args[1]         # known
    q = args[2]         # known
    gamma3 = args[3]    # known
    gamma4 = args[4]    # known
    dt = args[5]        # known
    dT_deriv = args[6]  # known
    r_t = dp_model_recomp_shared(k, t, q, gamma4, dt, r0)
    d_rt = forward_first_derivative(r_t, dt)
    tt = t[1:]
    rr = r_t[1:]
    dT_deriv_comp = (q/(4*pi*k*tt*rr))*np.exp(gamma3[1:]/tt)*(rr-2.0*tt*d_rt)
    diff = (dT_deriv-dT_deriv_comp)**2
    dsum = np.sum(diff)
    return dsum


def recomp_radius(fs, q, k, t, dT_known, r0_start=R0_START):
    dt = 1 / fs
    term = -(4 * pi * k) / q
    gamma2 = term * dT_known
    gamma3 = ei_inv_neg_array(gamma2) * t
    gamma4 = -gamma3
    dT_deriv = forward_first_derivative(dT_known, dt)
    args = (k, t, q, gamma3, gamma4, dt, dT_deriv)
    r0_start_in = np.asarray((r0_start,))
    res = minimize(_recomp, r0_start_in, args, method='Nelder-Mead', tol=1.0e-30)
    r0 = res.x[0]
    r_t = dp_model_recomp_shared(k, t, q, gamma4, dt, r0)
    return r0, r_t


def test_small_argument():
    from ei import ei
    from get_alpha import get_alpha_k_C
    from get_volumetric_heat_capacity import get_volumetric_heat_capacity
    from gen_time_vec import gen_time_vec

    q = 50
    k = 4
    r0 = 5.0e-3
    r1 = 50e-3
    theta_o = 1.0e-3
    theta_w = 0.49
    theta_m = 1 - theta_o - theta_w
    C = get_volumetric_heat_capacity(theta_m, theta_o, theta_w)
    fs = 120
    dt = 1 / fs
    T = 30
    t = gen_time_vec(fs, T)
    t = t[1:]  # cannot start with zero
    n = length(t)
    r = np.linspace(r0, r1, n)
    alpha = get_alpha_k_C(k, C)
    print('alpha = ', alpha)
    gamma4 = r**2 / (4 * alpha)

    gamma5 = np.log(np.sqrt(gamma4))
    gamma5_comp = np.log(r) - np.log(2*np.sqrt(alpha))
    gamma5_diff = gamma5 - gamma5_comp

    # take the average
    gamma6 = np.average(gamma5)
    N = length(gamma5)
    gamma6_comp = (1.0/N)*np.sum(np.log(r)) - np.log(2*np.sqrt(alpha))
    print('gamma6 = ', gamma6)
    print('gamma6_comp = ', gamma6_comp)
    print('diff = ', gamma6-gamma6_comp)

    # subtract from original
    gamma7 = gamma5 - gamma6
    gamma7_comp = np.log(r) - (1/N)*(np.sum(np.log(r)))
    gamma7_diff = gamma7 - gamma7_comp

    # # gamma7_comp_alt = np.log(r) - (1/N)*np.sum(np.log(r))
    # # gamma7_comp_alt = np.log(r) - np.sum(np.log(r**(1/N)))
    # gamma7_comp_alt = np.log(r/np.prod(r**(1/N)))
    # gamma7_diff_alt = gamma7 - gamma7_comp_alt
    #
    # gamma8 = np.exp(gamma7)
    # gamma8_comp = r/np.prod((r**(1/N)))
    #
    # gamma9 = gamma8*0.01
    # gamma9_diff = gamma9 - r
    #
    # # # check to see what this would be
    # # log_r_1overN = np.log(r**(1/N))
    #
    # arr = r**(1/N)
    # prod = np.prod(arr)
    # print('prod = ', prod)
    # print('N = ', N)
    #
    # av_log = np.mean(np.log(r))
    # print('av_log = ', av_log)
    #
    # block = False
    # plt.figure()
    # plt.plot(np.log(r))
    # plt.title('log')
    # plt.show(block)
    #
    # block = True
    # plt.figure()
    # plt.plot(gamma7)
    # plt.title('gamma7')
    # plt.show(block)

    # block = False
    # plt.figure()
    # plt.plot(gamma5)
    # plt.plot(gamma5_comp)
    # plt.title('Gamma5 comparison')
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(gamma5_diff)
    # plt.title('Gamma5 difference')
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(gamma7)
    # plt.plot(gamma7_comp)
    # plt.title('Gamma7 comparison')
    # plt.show(block)
    #
    # block = False
    # plt.figure()
    # plt.plot(gamma7_diff)
    # plt.title('Gamma7 difference')
    # plt.show(block)

    # block = True
    # plt.figure()
    # plt.plot(gamma7)
    # plt.plot(gamma7_comp_alt)
    # plt.title('Gamma7 alternate')
    # plt.show(block)

    # block = True
    # plt.figure()
    # plt.plot(gamma8)
    # plt.plot(gamma8_comp)
    # plt.title('Gamma8')
    # plt.show(block)

    # block = False
    # plt.figure()
    # plt.plot(gamma9)
    # plt.plot(r)
    # plt.title('Gamma9')
    # plt.show(block)
    #
    # block = True
    # plt.figure()
    # plt.plot(gamma9_diff)
    # plt.title('Gamma9 diff')
    # plt.show(block)


def test_full_inverse_heating():
    q = 50
    k = 4
    r0 = 5.7e-3
    r1 = 15e-3
    theta_o = 1.0e-3
    theta_w = 0.49
    theta_m = 1 - theta_o - theta_w
    C = get_volumetric_heat_capacity(theta_m, theta_o, theta_w)
    fs = 12
    dt = 1/fs
    T = 30
    t = gen_time_vec(fs, T)
    t = t[1:]   # cannot start with zero
    n = length(t)
    r = np.linspace(r0, r1, n)
    alpha = get_alpha_k_C(k, C)
    t_heating = 8
    th = t_heating

    # dT = compute_delta_T_dual_infinite_line_heating(q, k, alpha, r, t)  # obtain the forward model
    # r0_found, r_t = recomp_radius(fs, q, k, t, dT)
    # print('r0 = ', r0)
    # print('r0_found = ', r0_found)

    gamma4 = (r**2)/(4*alpha)

    gamma5 = np.sqrt(gamma4)
    gamma5_comp = r/(2*np.sqrt(alpha))

    gamma6 = np.log(gamma5)
    gamma7 = forward_first_derivative(gamma6, dt)
    gamma7_comp = forward_first_derivative(np.log(r), dt)

    r0_in = r0
    gamma8 = np.exp(inverse_first_derivative(gamma7, dt, np.log(r0_in)))
    gamma8 = gamma8 - gamma8[0]
    print('delta_r = ', r[2]-r[1])
    print('delta_r = ', gamma8[2] - gamma8[1])

    plt.figure()
    plt.plot(gamma8)
    plt.plot(r)
    plt.show()


def test_matrix_inverse():
    """
    Test the matrix inverse
    :return:
    """
    A = np.zeros((3, 3))
    N = 3

    A[0, 0] = 1
    A[0, 1] = 0
    A[0, 2] = 0

    A[1, 0] = 0
    A[1, 1] = 1
    A[1, 2] = 0

    A[2, 0] = 0
    A[2, 1] = 0
    A[2, 2] = 1
    print(np.linalg.det(A))


from scipy.optimize import minimize
def test_full_inverse_heating_determine_r():
    q = 50
    k = 4
    r0 = 5.7e-3
    r1 = 15e-3
    theta_o = 1.0e-3
    theta_w = 0.49
    theta_m = 1 - theta_o - theta_w
    C = get_volumetric_heat_capacity(theta_m, theta_o, theta_w)
    fs = 12
    dt = 1/fs
    T = 30
    t = gen_time_vec(fs, T)
    t = t[1:]           # cannot start with zero
    n = length(t)
    r = np.linspace(r0, r1, n)
    alpha = get_alpha_k_C(k, C)
    t_heating = 8
    th = t_heating

    # this is the known curve
    print('Running for different starting r0 values')
    dT = compute_delta_T_dual_infinite_line_heating(q, k, alpha, r, t)
    r0v = [4e-3, 5e-3, 5.7e-3, 6e-3]
    dT_inv_vec = []
    for r0_elem in r0v:
        print('r0_elem = ', r0_elem)
        r_t_det = dp_model_recomp_heating(fs, q, k, t, dT, r0_elem, get_gamma4=False)
        dT_inv = compute_delta_T_dual_infinite_line_heating(q, k, alpha, r_t_det, t)
        dT_inv_vec.append(dT_inv)

    print('Starting optimization to find the known r')
    rstart = 1.0e-3

    # def _f_obtain_r_from_curve(x, *args):
    #     # UNKNOWN
    #     r0_in = x[0]
    #     # KNOWN
    #     q_in = args[0]
    #     k_in = args[1]
    #     alpha_in = args[2]
    #     t_in = args[3]
    #     fs_in = args[4]
    #     known_in = args[5]
    #     r_t_det = dp_model_recomp_heating(fs, q_in, k_in, t_in, known_in, r0_in, get_gamma4=False)
    #     dT_inv = compute_delta_T_dual_infinite_line_heating(q_in, k_in, alpha_in, r_t_det, t_in)
    #     diff = (known_in - dT_inv)**2
    #     dsum = np.sum(diff)
    #     return dsum
    # # use optimization for the starting r
    # xstart = (rstart,)
    # args = (q, k, alpha, t, fs, dT)
    # res = minimize(_f_obtain_r_from_curve, xstart, args, method='Nelder-Mead')
    # rout = res.x[0]
    # print('rout = ', rout)
    # print('r0 = ', r0)
    #
    block = True
    plt.figure()
    plt.plot(t, dT, linewidth=8)
    for dT_i in dT_inv_vec:
        plt.plot(t, dT_i)
    plt.show(block)


def main():
    # test_full_inverse_heating()
    # test_matrix_inverse()
    test_full_inverse_heating_determine_r()


if __name__ == '__main__':
    main()



