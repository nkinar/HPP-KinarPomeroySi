import numpy as np
from calibration_processing import calibration_processing
from constant_labels import create_label
from ticklabels_off import ticklabels_off_x
from comparisons import compute_rmse, compute_mb, compute_percentage_diff
from place_rmse_mb_on_plot import place_rmse_mb_pd_on_plot
from set_xlim_linear import set_xlim_linear
from EasyDump import EasyDump
from obtain_load_r0 import SaveLoadR0
from constants import *

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('ggplot')
mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams.update({'font.size': FONT_SIZE})

"""
File to obtain calibration figures for the paper
"""


def calibration_plot_table(path, q_example, show):
    """
    Create the calibration plot and table showing how different q values are affected
    :param path:            as the path to the calibration data
    :param q_list:          as the list of q values over which to iterate
    :param q_example:       as the q for which to plot an example
    :param show:            True to show the figure
    :return:
    """
    q = q_example
    additional_text = NO_TEXT

    #########################################################
    t_cold = T_COLD_BEGIN
    t_heat = T_HEAT_CAL             # seconds heating (s)
    t_total = TIME_TOTAL_CAL        # 3 minutes
    fs = FS_SAMPLE                  # sampling rate (Hz)
    k_assumed = K_WATER             # thermal conductivity of water (W m^-1 K^-1)
    rho = RHO_WATER_AGAR            # density of agar gel water (kg m^-3)
    c = HEAT_CAPACITY_WATER         # heat capacity of water (J kg^-1 K^-1)
    filt = 'both'                   # filter both SP and DP
    downsample = 'dp'               # downsample DP
    fs_down = DOWNSAMPLE_DP         # 1 second downsample
    filt_cutoff = DP_LOWPASS_FILTER_CUT  # lowpass filter cutoff frequency
    use_assumed_k = False
    use_assumed_q = False
    use_step_detector_q = True
    #########################################################

    # iterate over the variables to create the first plot
    use_assumed_k_range = [True, False]
    q_range = Q_VAL_CAL_LIST
    r0_assumed_k = []       # r0 with an assumed k
    r0_find_k = []          # r0 wtih a found k using the single probe
    found_k = []            # k determined using the single probe
    found_q = []            # q determined by averaging over the plateau

    sp_measured = []
    dp_measured = []

    sp_model_assumed_k = []
    sp_model_determined_k = []

    dp_model_assumed_k = []
    dp_model_determined_k = []

    #############################
    # STATISTICS
    #############################
    for use_assumed_k in use_assumed_k_range:
        print('use_assumed_k: ', use_assumed_k)
        if use_assumed_k:
            r0_vec = r0_assumed_k
            sp_model = sp_model_assumed_k
            dp_model = dp_model_assumed_k
        else:
            r0_vec = r0_find_k
            sp_model = sp_model_determined_k
            dp_model = dp_model_determined_k
        for q in q_range:
            print('q = ', q)
            original, delta_T1, delta_T2, t1, t2, delta_T1_trim, delta_T2_trim, t1_trim, t2_trim, \
            delta_T1_trim_heating, t1_trim_heating, qav, step_q_tfirst, step_q_tsecond, \
            ypk, k, alpha, rfound, delta_T_single_synth, deltaT_dual_synth, t, rmse_single, mb_single, pd_single, \
            rmse_dual, mb_dual, pd_dual, rmse_q_calc, mb_q_calc, pd_q_calc, t2_trim_heating, delta_T2_trim_heating = \
            calibration_processing(path, additional_text, q, use_assumed_q, use_step_detector_q, t_cold, t_heat,
                                       t_total,
                                       fs, k_assumed,
                                       use_assumed_k, rho, c, downsample, fs_down, filt, filt_cutoff)
            r0_vec.append(rfound)                       # find the r value
            sp_model.append(delta_T_single_synth)       # find the single probe model
            dp_model.append(deltaT_dual_synth)          # find the dual probe model
            if not use_assumed_k:
                found_k.append(k)                       # k determined from the single probe
                # the following are the same when use_assumed_k==True, but we only want to calculate once
                # to do comparisons with the actual data
                found_q.append(qav)
                sp_measured.append(delta_T1_trim_heating[1:])
                dp_measured.append(delta_T2_trim[1:])

    #############################
    # STATISTICS
    #############################
    print('Writing statistics to file')
    ed = EasyDump(TABLE_PATH + CAL_PLOT_TABLE_FILENAME + TABLE_EXT, 2)
    ed.open()

    sd_r0_assumed_k = np.std(np.asarray(r0_assumed_k)/1.0e-3)  # in mm
    sd_r0_determined_k = np.std(np.asarray(r0_find_k)/1.0e-3)  # in mm
    ed.write(sd_r0_assumed_k, 'sd_r0_assumed_k')
    ed.write(sd_r0_determined_k, 'sd_r0_determined_k')

    rmse_k = compute_rmse(k_assumed, found_k)
    mb_k = compute_mb(k_assumed, found_k)
    pd_k = compute_percentage_diff(k_assumed, found_k)
    ed.write(rmse_k, 'rmse_k')
    ed.write(mb_k, 'mb_k')
    ed.write(pd_k, 'pd_k')

    rmse_q = compute_rmse(q_range, found_q)
    mb_q = compute_mb(q_range, found_q)
    pd_q = compute_percentage_diff(q_range, found_q)
    ed.write(rmse_q, 'rmse_q')
    ed.write(mb_q, 'mb_q')
    ed.write(pd_q, 'pd_q')

    sp_measured = np.vstack(sp_measured)
    dp_measured = np.vstack(dp_measured)
    sp_model_assumed_k = np.vstack(sp_model_assumed_k)
    sp_model_determined_k = np.vstack(sp_model_determined_k)
    dp_model_assumed_k = np.vstack(dp_model_assumed_k)
    dp_model_determined_k = np.vstack(dp_model_determined_k)

    rmse_sp_assumed_k = compute_rmse(sp_measured, sp_model_assumed_k)
    rmse_sp_determined_k = compute_rmse(sp_measured, sp_model_determined_k)
    rmse_dp_assumed_k = compute_rmse(dp_measured, dp_model_assumed_k)
    rmse_dp_determined_k = compute_rmse(dp_measured, dp_model_determined_k)

    mb_sp_assumed_k = compute_mb(sp_measured, sp_model_assumed_k)
    mb_sp_determined_k = compute_mb(sp_measured, sp_model_determined_k)
    mb_dp_assumed_k = compute_mb(dp_measured, dp_model_assumed_k)
    mb_dp_determined_k = compute_mb(dp_measured, dp_model_determined_k)

    pd_sp_assumed_k = compute_percentage_diff(sp_measured, sp_model_assumed_k)
    pd_sp_determined_k = compute_percentage_diff(sp_measured, sp_model_determined_k)
    pd_dp_assumed_k = compute_percentage_diff(dp_measured, dp_model_assumed_k)
    pd_dp_determined_k = compute_percentage_diff(dp_measured, dp_model_determined_k)

    ed.write(rmse_sp_assumed_k, 'rmse_sp_assumed_k')
    ed.write(mb_sp_assumed_k, 'mb_sp_assumed_k')
    ed.write(pd_sp_assumed_k, 'pd_sp_assumed_k')

    ed.write(rmse_sp_determined_k, 'rmse_sp_determined_k')
    ed.write(mb_sp_determined_k, 'mb_sp_determined_k')
    ed.write(pd_sp_determined_k, 'pd_sp_determined_k')

    ed.write(rmse_dp_assumed_k, 'rmse_dp_assumed_k')
    ed.write(mb_dp_assumed_k, 'mb_dp_assumed_k')
    ed.write(pd_dp_assumed_k, 'pd_dp_assumed_k')

    ed.write(rmse_dp_determined_k, 'rmse_dp_determined_k')
    ed.write(mb_dp_determined_k, 'mb_dp_determined_k')
    ed.write(pd_dp_determined_k, 'pd_dp_determined_k')

    # close the file
    ed.close()
    print('DONE writing statistics to file')

    #############################
    # EXAMPLE at a fixed q value
    #############################
    original, delta_T1, delta_T2, t1, t2, delta_T1_trim, delta_T2_trim, t1_trim, t2_trim, \
    delta_T1_trim_heating, t1_trim_heating, qav, step_q_tfirst, step_q_tsecond, \
    ypk, k, alpha, rfound, delta_T_single_synth, deltaT_dual_synth, t, rmse_single, mb_single, pd_single, \
    rmse_dual, mb_dual, pd_dual, rmse_q_calc, mb_q_calc, pd_q_calc, t2_trim_heating, delta_T2_trim_heating = \
        calibration_processing(path, additional_text, q_example, use_assumed_q, use_step_detector_q, t_cold,
                               t_heat, t_total, fs, k_assumed, use_assumed_k, rho, c,
                               downsample, fs_down, filt, filt_cutoff)
    num, Vout, dV, I, Rw, T1, T2, dac_code, qprime = original  # unpack the original data

    r0_assumed_k_mm = np.asarray(r0_assumed_k)/1.0e-3
    r0_find_k_mm = np.asarray(r0_find_k)/1.0e-3

    #################################
    # PLOTS
    #################################

    fig = plt.figure(num=None, figsize=SQUARE_FIGSIZE)
    nx = 6      # number of ticks on the x axis for the time

    # NOTE that rmse_single, mb_single, rmse_dual, mb_dual are to be shown on the plots using the text

    # plot of the r0
    ax = fig.add_subplot(2, 2, 1)
    ax.plot(q_range, r0_assumed_k_mm, 'bs', label='Assumed $k$')
    ax.plot(q_range, r0_find_k_mm, 'g^', label='Determined $k$')
    ax.set_ylim([8, 11])
    set_xlim_linear(q_range[0], q_range[-1], nx)
    ax.set_xlim([q_range[0]-3, q_range[-1]+3])
    ax.legend()
    ax.set_xlabel(create_label('$q$', 'W m^-1'))
    ax.set_ylabel(create_label('$r_0$', 'mm'))
    ax.set_title('(a)', loc='left')

    # EXAMPLE plot
    use_pd = True
    is_scientific = True

    # plot of single probe example (SP)
    ax = fig.add_subplot(2, 2, 2)
    ax.plot(t1_trim, delta_T1_trim, color=MEASURED_COLOR, label='Measured')
    ax.plot(t1_trim_heating[1:], delta_T_single_synth, color=MODELLED_COLOR, label='Modelled')
    place_rmse_mb_pd_on_plot(rmse_single, mb_single, pd_single, 'K', 1, 'center right', is_scientific, use_pd)
    set_xlim_linear(0, t_total, nx)
    ax.set_xlim([-3, t_total+3])
    ax.legend()
    ax.set_xlabel('Time (s)')
    ticklabels_off_x()  # do not show tick labels for this plot
    ax.set_ylabel(create_label('$\Delta T$', 'K'))
    ax.set_title('(b)', loc='left')

    # plot of q over the time of heating
    ax = fig.add_subplot(2, 2, 3)
    ax.plot(t, qprime, color=MEASURED_COLOR, label='Measured')
    ax.axvline(x=step_q_tfirst, linestyle=':', color=MODELLED_COLOR, label='Step Detector')
    ax.axvline(x=step_q_tsecond, linestyle=':', color=MODELLED_COLOR)
    ax.legend()
    place_rmse_mb_pd_on_plot(rmse_q_calc, mb_q_calc, pd_q_calc, 'W \hspace{0.1} m^{-1}', 1, 'upper center',
                             is_scientific,
                             use_pd)
    ax.set_xlim([0, t_heat+2])
    ax.set_ylim([-1, q+10])
    ax.set_xlabel('Time(s)')
    ax.set_ylabel(create_label('$q$', 'W \hspace{0.1} m^-1'))
    ax.set_title('(c)', loc='left')

    # plot of dual probe example (DP)
    ax = fig.add_subplot(2, 2, 4)
    ax.plot(t2_trim, delta_T2_trim, color=MEASURED_COLOR, label='Measured')
    ax.plot(t2_trim[1:], deltaT_dual_synth, color=MODELLED_COLOR, label='Modelled')
    place_rmse_mb_pd_on_plot(rmse_dual, mb_dual, pd_dual, 'K', 1, 'lower right', is_scientific,
                             use_pd)
    set_xlim_linear(0, t_total, nx)
    ax.set_xlim([-3, t_total+3])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(create_label('$\Delta T$', 'K'))
    ax.set_title('(d)', loc='left')

    # save the plot out to file
    plt.tight_layout()
    plt.savefig(FIG_PATH + CAL_PLOT_FIGURE_FILENAME + PLOT_EXTENSION)

    # save out the r0
    sr = SaveLoadR0()
    sr.save_r0(CONST_PATH + CAL_R0_FILENAME + SAVEZ_EXT, Q_VAL_CAL_LIST, r0_find_k_mm)

    if show:   # show the plot for the purposes of debugging
        plt.show()
    # DONE


def run_calibration_plot():
    path = '../hpp-data/hpp-formatted-data/calibration/July-12-2017/'
    q_list = []
    q_example = 50
    show = False
    calibration_plot_table(path, q_example, show)


def main():
    print('Running calibration plot')
    run_calibration_plot()
    print('DONE running calibration plot')


if __name__ == '__main__':
    main()

