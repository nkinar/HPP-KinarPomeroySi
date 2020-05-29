import numpy as np
from calibration_processing import initial_compute_hpp
from comparisons import compute_rmse, compute_mb
from inverse_model_nominal import obtain_sp_vars_from_curve_min, compute_single_probe_forward_model_min
from curve_time_trim import curve_time_trim
from inverse_model_nominal import obtain_k_H_from_dual_probe, obtain_sp_vars_from_curve
from sp_model_signal import sp_model_late_time_inv, sp_model_signal_inv
from comparisons import compute_rmse, compute_mb, compute_percentage_diff
from EasyDump import EasyDump
from constant_labels import create_label
from set_xlim_linear import set_xlim_linear
from constants import *

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('ggplot')
mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams.update({'font.size': FONT_SIZE})


def test_calibration_k(path, start_string, additional_text, downsample, filt, filt_cutoff, fs,
                         fs_down, q, t_cold, t_heat, t_total, k_assumed, Hstart):
    """
    Obtain the k from a series of files and test to see how the k can be determined
    :param path:
    :param start_string:
    :param additional_text:
    :param downsample:
    :param filt:
    :param filt_cutoff:
    :param fs:
    :param fs_down:
    :param q:
    :param t_cold:
    :param t_heat:
    :param t_total:
    :param k_assumed:
    :param Hstart:
    :return:
    """
    use_assumed_q = False
    use_step_detector_q = True
    I, Rw, T1, T2, Vout, dV, dac_code, delta_T1, delta_T1_trim, delta_T1_trim_heating, delta_T2, \
    delta_T2_trim, num, qav, qprime, step_q_tfirst, step_q_tsecond, t1, t1_trim, t1_trim_heating, t2, t2_trim, \
    ypk, t, rmse_q_calc, mb_q_calc, pd_q_calc, t2_trim_heating, delta_T2_trim_heating = \
        initial_compute_hpp(start_string, additional_text,
                            downsample, filt, filt_cutoff, fs, fs_down, path, q, t_cold, t_heat, t_total,
                            use_assumed_q,
                            use_step_detector_q)

    # use the average q determined from the curve
    qq = qav
    # starting k for inverse models
    kstart = k_assumed

    # CURVE-FITTING HEATING DP
    t2_trim_heating1, delta_T2_trim_heating1 = curve_time_trim(t2_trim_heating, delta_T2_trim_heating)
    t_heating = t_heat
    full = False
    kout_dp_heating, _hout = obtain_k_H_from_dual_probe(qq, t2_trim_heating1, t_heating, delta_T2_trim_heating1, kstart,
                                                        Hstart,
                                                        full)

    # CURVE-FITTING HEATING AND COOLING
    t2_trim1, delta_T2_trim1 = curve_time_trim(t2_trim, delta_T2_trim)
    full = True
    kout_dp_heating_cooling, _hout = obtain_k_H_from_dual_probe(q, t2_trim1, t_heating, delta_T2_trim1, kstart, Hstart,
                                                                full)

    # LINEAR LATE-TIME MODEL
    t1_trim_heating1, delta_T1_trim_heating1 = curve_time_trim(t1_trim_heating, delta_T1_trim_heating)
    dt = t1_trim_heating[1] - t1_trim_heating[0]
    tlinear, kout_sp_late_time = sp_model_late_time_inv(delta_T1_trim_heating1, t1_trim_heating1, dt, qav)
    # NOTE that tlinear is with respect to the cut curve

    # LINEAR LATE-TIME MODEL (ALL DATA)
    entire = True
    _tlinear, kout_sp_all = sp_model_late_time_inv(delta_T1_trim_heating1, t1_trim_heating1, dt, qav, entire)

    # SIGNAL PROCESSING TO OBTAIN K
    kout_signal, bdet_signal = sp_model_signal_inv(delta_T1_trim_heating1, t1_trim_heating1, dt, qav)

    # CURVE-FITTING AFTER SIGNAL PROCESSING
    kout_signal_curve_fit, _bd, _cd, _dd = obtain_sp_vars_from_curve(q, t1_trim_heating1, delta_T1_trim_heating1,
                                                                     kout_signal)
    """
    OUTPUTS: 
    kout_dp_heating             = curve-fit to dual probe heating
    kout_dp_heating_cooling     = curve-fit to dual probe heating and cooling
    kout_sp_late_time           = single probe late-time model
    tlinear                     = time detected when the curve becomes linear
    kout_sp_all                 = single probe late-time model (all heating data)
    kout_signal                 = signal processing to obtain k 
    bdet_signal                 = b that is detected from the model
    kout_signal_curve_fit       = signal processing and then fitting to the full SP model
    """
    return kout_dp_heating, kout_dp_heating_cooling, kout_sp_late_time, tlinear, kout_sp_all, \
            kout_signal, bdet_signal, kout_signal_curve_fit


def run_test_calibration_k(show_plot):
    path = '../hpp-data/hpp-formatted-data/calibration/July-12-2017/'
    dec_places = 2                              # number of decimal places to have result
    t_cold = T_COLD_BEGIN                       # time at beginning which is cold
    t_heat = T_HEAT_CAL                         # seconds heating (s)
    t_total = TIME_TOTAL_CAL                    # 3 minutes for entire test
    fs = FS_SAMPLE                              # sampling rate (Hz)
    k_assumed = K_WATER                         # thermal conductivity of water (W m^-1 K^-1)
    rho = RHO_WATER_AGAR                        # density of agar gel water (kg m^-3)
    c = HEAT_CAPACITY_WATER                     # heat capacity of water (J kg^-1 K^-1)
    filt = 'both'                               # filter both SP and DP
    downsample = 'dp'                           # downsample DP
    fs_down = DOWNSAMPLE_DP                     # 1 second downsample
    filt_cutoff = DP_LOWPASS_FILTER_CUT         # lowpass filter cutoff frequency
    start_string = CAL_BEGIN                    # start calibration string
    additional_text = NO_TEXT                   # additional text
    Hstart = 80                                 # starting value for H
    qlist = Q_VAL_CAL_LIST                      # list of q values during calibration

    kout_dp_heating_vec = []
    kout_dp_heating_cooling_vec = []
    kout_sp_late_time_vec = []
    kout_sp_all_vec = []
    kout_signal_vec = []
    kout_signal_curve_fit_vec = []

    for q in qlist:
        print('q = ', q)
        kout_dp_heating, kout_dp_heating_cooling, kout_sp_late_time, _tlinear, kout_sp_all, \
        kout_signal, bdet_signal, kout_signal_curve_fit = \
        test_calibration_k(path, start_string, additional_text, downsample, filt, filt_cutoff, fs,
                           fs_down, q, t_cold, t_heat, t_total, k_assumed, Hstart)

        kout_dp_heating_vec.append(kout_dp_heating)
        kout_dp_heating_cooling_vec.append(kout_dp_heating_cooling)
        kout_sp_late_time_vec.append(kout_sp_late_time)
        kout_sp_all_vec.append(kout_sp_all)
        kout_signal_vec.append(kout_signal)
        kout_signal_curve_fit_vec.append(kout_signal_curve_fit)
    # DONE

    qlist = np.asarray(qlist)
    kout_dp_heating_vec = np.asarray(kout_dp_heating_vec)
    kout_dp_heating_cooling_vec = np.asarray(kout_dp_heating_cooling_vec)
    kout_sp_late_time_vec = np.asarray(kout_sp_late_time_vec)
    kout_sp_all_vec = np.asarray(kout_sp_all_vec)
    kout_signal_vec = np.asarray(kout_signal_vec)
    kout_signal_curve_fit_vec = np.asarray(kout_signal_curve_fit_vec)

    ed = EasyDump(TABLE_PATH + K_PLOT_TABLE_FILENAME + TABLE_EXT, 2)
    ed.open()

    # COMPUTE RMSE AND MB FOR ALL EXPERIMENTS
    rmse_kout_dp_heating = compute_rmse(k_assumed, kout_dp_heating_vec)
    mb_kout_dp_heating = compute_mb(k_assumed, kout_dp_heating_vec)
    pd_kout_dp_heating = compute_percentage_diff(k_assumed, kout_dp_heating_vec)
    ed.write(rmse_kout_dp_heating, 'rmse_kout_dp_heating')
    ed.write(mb_kout_dp_heating, 'mb_kout_dp_heating')
    ed.write(pd_kout_dp_heating, 'pd_kout_dp_heating')

    rmse_kout_dp_heating_cooling = compute_rmse(k_assumed, kout_dp_heating_cooling_vec)
    mb_kout_dp_heating_cooling = compute_mb(k_assumed, kout_dp_heating_cooling_vec)
    pd_kout_dp_heating_cooling = compute_percentage_diff(k_assumed, kout_dp_heating_cooling_vec)
    ed.write(rmse_kout_dp_heating_cooling, 'rmse_kout_dp_heating_cooling')
    ed.write(mb_kout_dp_heating_cooling, 'mb_kout_dp_heating_cooling')
    ed.write(pd_kout_dp_heating_cooling, 'pd_kout_dp_heating_cooling')

    rmse_kout_sp_late_time = compute_rmse(k_assumed, kout_sp_late_time_vec)
    mb_kout_sp_late_time = compute_mb(k_assumed, kout_sp_late_time_vec)
    pd_kout_sp_late_time = compute_percentage_diff(k_assumed, kout_sp_late_time_vec)
    ed.write(rmse_kout_sp_late_time, 'rmse_kout_sp_late_time')
    ed.write(mb_kout_sp_late_time, 'mb_kout_sp_late_time')
    ed.write(pd_kout_sp_late_time, 'pd_kout_sp_late_time')

    rmse_kout_sp_all = compute_rmse(k_assumed, kout_sp_all_vec)
    mb_kout_sp_all = compute_mb(k_assumed, kout_sp_all_vec)
    pd_kout_sp_all = compute_percentage_diff(k_assumed, kout_sp_all_vec)
    ed.write(rmse_kout_sp_all, 'rmse_kout_sp_all')
    ed.write(mb_kout_sp_all, 'mb_kout_sp_all')
    ed.write(pd_kout_sp_all, 'pd_kout_sp_all')

    rmse_kout_signal = compute_rmse(k_assumed, kout_signal_vec)
    mb_kout_signal = compute_mb(k_assumed, kout_signal_vec)
    pd_kout_signal = compute_percentage_diff(k_assumed, kout_signal_vec)
    ed.write(rmse_kout_signal, 'rmse_kout_signal')
    ed.write(mb_kout_signal, 'mb_kout_signal')
    ed.write(pd_kout_signal, 'pd_kout_signal')

    rmse_kout_signal_curve_fit = compute_rmse(k_assumed, kout_signal_curve_fit_vec)
    mb_kout_signal_curve_fit = compute_mb(k_assumed, kout_signal_curve_fit_vec)
    pd_kout_signal_curve_fit = compute_percentage_diff(k_assumed, kout_signal_curve_fit_vec)
    ed.write(rmse_kout_signal_curve_fit, 'rmse_kout_signal_curve_fit')
    ed.write(mb_kout_signal_curve_fit, 'mb_kout_signal_curve_fit')
    ed.write(pd_kout_signal_curve_fit, 'pd_kout_signal_curve_fit')

    # close the file
    ed.close()

    # need to trim outliers
    elem = kout_dp_heating_vec < 100
    qq = qlist[elem]
    kout_dp_heating_vec_trim = kout_dp_heating_vec[elem]

    elem = kout_dp_heating_vec_trim > 0.01
    qq = qq[elem]
    kout_dp_heating_vec_trim = kout_dp_heating_vec_trim[elem]

    fig = plt.figure(num=None)
    nx = 6
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(qq, kout_dp_heating_vec_trim, "o",  label='Heating DP', color='orange')
    ax.plot(qlist, kout_dp_heating_cooling_vec, "v", label='Heating and Cooling DP')
    ax.plot(qlist, kout_sp_late_time_vec, "s", label='Late-Time SP')
    ax.plot(qlist, kout_sp_all_vec, "^", label='Late-Time SP (all data)')
    ax.plot(qlist, kout_signal_vec, "*", label='Signal Processing SP', color="yellowgreen")
    ax.plot(qlist, kout_signal_curve_fit_vec, "h", label='Signal Processing SP with Curve-Fitting')
    ax.set_ylim([-0.1, 7])
    set_xlim_linear(qlist[0], qlist[-1], nx)
    ax.set_xlim([qlist[0]-3, qlist[-1]+3])
    ax.set_xlabel(create_label('$q$', 'W \hspace{0.1} m^-1'))
    ax.set_ylabel(create_label('$k$', 'W \hspace{0.1} m^-1 \hspace{0.1} K^-1'))
    ax.axhline(y=k_assumed, color="silver", linestyle='dashed', label='Assumed Value')
    ax.legend(loc='upper left')

    plt.tight_layout()
    plt.savefig(FIG_PATH + K_PLOT_FILENAME + PLOT_EXTENSION)

    if show_plot:   # show the plot for purposes of debugging
        plt.show()


def main():
    show_plot = False
    run_test_calibration_k(show_plot)


if __name__ == '__main__':
    main()
