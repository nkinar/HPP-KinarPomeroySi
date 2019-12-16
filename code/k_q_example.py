import numpy as np
from calibration_processing import initial_compute_hpp
from curve_time_trim import curve_time_trim
from sp_model_signal import sp_model_late_time_inv, sp_model_signal_inv, sp_model_nominal_get_b_c_d
from constant_labels import create_label
from sp_model_signal import sp_model_late_time, dT_sp_model_signal
from comparisons import compute_percentage_diff, compute_mb, compute_rmse
from place_rmse_mb_on_plot import place_rmse_mb_pd_on_plot
from constants import *
from numpy import pi

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('ggplot')
mpl.rcParams["mathtext.fontset"] = "stix"


def run_k_q_example(start_string, additional_text, downsample, filt, filt_cutoff, fs, fs_down,
                    path, q, t_cold, t_heat, t_total, show, fn):
    use_assumed_q = False
    use_step_detector_q = True
    I, Rw, T1, T2, Vout, dV, dac_code, delta_T1, delta_T1_trim, delta_T1_trim_heating, \
    delta_T2, delta_T2_trim, num, qav, qprime, step_q_tfirst, step_q_tsecond, t1, t1_trim, \
    t1_trim_heating, t2, t2_trim, ypk, t, rmse_q_calc, mb_q_calc, pd_q_calc, t2_trim_heating, \
    delta_T2_trim_heating = \
                            initial_compute_hpp(start_string, additional_text,
                            downsample, filt, filt_cutoff, fs, fs_down, path, q, t_cold, t_heat, t_total,
                            use_assumed_q,
                            use_step_detector_q)
    t1_trim_heating1, delta_T1_trim_heating1 = curve_time_trim(t1_trim_heating, delta_T1_trim_heating)
    dt = t1_trim_heating1[1] - t1_trim_heating1[0]

    # late-time model
    entire = False
    return_all = True
    tlinear, kdet, bdet, t_tlinear0, dT_tlinear0 = sp_model_late_time_inv(delta_T1_trim_heating1,
                                                                    t1_trim_heating1, dt, qav, entire, return_all)
    tlog = np.log(t1_trim_heating1)
    dT_sp_synth = sp_model_late_time(qav, kdet, t_tlinear0, bdet)
    rmse_sp = compute_rmse(dT_tlinear0, dT_sp_synth)
    mb_sp = compute_mb(dT_tlinear0, dT_sp_synth)
    pd_sp = compute_percentage_diff(dT_tlinear0, dT_sp_synth)

    # signal processing model
    kdet_sig, bdet_sig, tth, dT2_h = sp_model_signal_inv(delta_T1_trim_heating1, t1_trim_heating1, dt, q,
                                                         output_model_array=True)
    bout, cout, dout = sp_model_nominal_get_b_c_d(delta_T1_trim_heating1, t1_trim_heating1, q, kdet_sig)
    print('SP output:')
    print('kdet = ', kdet)
    print('bout = ', bout)
    print('cout = ', cout)
    print('dout = ', dout)

    # synthetic signal processing for comparison
    dT_model_sig_synth = dT_sp_model_signal(qav, kdet_sig, tth, bdet_sig)

    rmse_sp_sig = compute_rmse(dT2_h, dT_model_sig_synth)
    mb_sp_sig = compute_mb(dT2_h, dT_model_sig_synth)
    pd_sp_sig = compute_percentage_diff(dT2_h, dT_model_sig_synth)

    fig = plt.figure(num=None, figsize=HSQUARE_FIGSIZE)
    is_scientific = True
    use_pd = True

    # (a)
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(tlog, delta_T1_trim_heating1, color=MEASURED_COLOR, label='Measured')
    ax.plot(np.log(t_tlinear0), dT_sp_synth, color=MODELLED_COLOR, label='Modelled')
    ax.axvline(x=np.log(tlinear), linestyle=':', color=MODELLED_COLOR, label='Identified Linear Section')
    place_rmse_mb_pd_on_plot(rmse_sp, mb_sp, pd_sp, 'K', 1, 'center left',
                             is_scientific,
                             use_pd)
    ax.legend()
    ax.set_xlabel(r'$\mathrm{log}(t)$')
    ax.set_ylabel(create_label('$\Delta T$', 'K'))
    ax.set_title('(a)', loc='left')

    # (b)
    ax = fig.add_subplot(1, 2, 2)
    ax.plot(tth, dT2_h, color=MEASURED_COLOR, label='signal')
    ax.plot(tth, dT_model_sig_synth, color=MODELLED_COLOR, label='model')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(create_label('$\Delta \Gamma_{5} \hspace{0.3} (\hspace{0.3}t\hspace{0.3})$', ''))
    place_rmse_mb_pd_on_plot(rmse_sp_sig, mb_sp_sig, pd_sp_sig, 'K', 1, 'lower right',
                             is_scientific,
                             use_pd)
    ax.set_title('(b)', loc='left')

    plt.tight_layout()
    plt.savefig(FIG_PATH + fn + PLOT_EXTENSION)

    if show:
        plt.show()


def run_k_example(show):
    path = '../hpp-data/hpp-formatted-data/calibration/July-12-2017/'
    q = 50

    t_cold = T_COLD_BEGIN                       # time at beginning that is cold
    t_heat = T_HEAT_CAL                         # seconds heating (s)
    t_total = TIME_TOTAL_CAL                    # 3 minutes for entire test
    fs = FS_SAMPLE                              # sampling rate (Hz)
    filt = 'both'                               # filter both SP and DP
    downsample = 'dp'                           # downsample DP
    fs_down = DOWNSAMPLE_DP                     # 1 second downsample
    filt_cutoff = DP_LOWPASS_FILTER_CUT         # lowpass filter cutoff frequency
    start_string = CAL_BEGIN                    # start calibration string
    additional_text = NO_TEXT                   # additional text
    run_k_q_example(start_string, additional_text, downsample, filt, filt_cutoff, fs, fs_down,
                    path, q, t_cold, t_heat, t_total, show, K_PLOT_FILENAME_EXAMPLE)


def run_k_example_sand(show):
    path = '../hpp-data/hpp-formatted-data/sand/July-20-2017/'
    q = 45

    t_cold = T_COLD_BEGIN                       # time at beginning that is cold
    t_heat = 8                                  # seconds heating (s)
    t_total = 3*60                              # 3 minutes for entire test
    fs = FS_SAMPLE                              # sampling rate (Hz)
    filt = 'both'                               # filter both SP and DP
    downsample = 'dp'                           # downsample DP
    fs_down = DOWNSAMPLE_DP                     # 1 second downsample
    filt_cutoff = DP_LOWPASS_FILTER_CUT         # lowpass filter cutoff frequency
    start_string = 'sand-'                      # start string to load the data
    additional_text = '-rep1'                   # additional text
    run_k_q_example(start_string, additional_text, downsample, filt, filt_cutoff, fs, fs_down,
                    path, q, t_cold, t_heat, t_total, show, K_PLOT_FILENAME_EXAMPLE_SAND)


def main():
    show_plot = False
    run_k_example(show_plot)
    run_k_example_sand(show_plot)


if __name__ == '__main__':
    main()
