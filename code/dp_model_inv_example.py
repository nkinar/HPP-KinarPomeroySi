import numpy as np
from forward_model_nominal import compute_delta_T_dual_infinite_line_heating, compute_delta_T_dual_heating_cooling
from get_alpha import get_alpha_k_C
from get_volumetric_heat_capacity import get_volumetric_heat_capacity
from gen_time_vec import gen_time_vec
from get_size import length
from constant_labels import create_label
from dp_model_recomp import dp_model_recomp_heating, dp_model_recomp_heating_cooling
from ticklabels_off import ticklabels_off_x
from brownian import brownian_noise_norm
from constants import *

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('ggplot')
mpl.rcParams["mathtext.fontset"] = "stix"


def run_plot_output(alpha, fs, k, q, r0, run_full, t_heat, t, r):
    """
    Function required to generate the plot outputs for the demonstration
    :param alpha:
    :param fs:
    :param k:
    :param q:
    :param r0:
    :param run_full:
    :param t_heat:
    :param t:
    :param r:
    :return:
    """
    # Note that Ei(x=0) is not defined, so the time vector must be trimmed
    # The forward model is not defined at t = 0.
    tt = np.copy(t)
    tt = tt[1:]
    dt = 1 / fs
    n = length(tt)
    if run_full:  # run the complete algorithm for heating and cooling
        dT_dual = compute_delta_T_dual_heating_cooling(q, k, alpha, r, tt, t_heat, split=False)
        dT_dual_fixed = compute_delta_T_dual_heating_cooling(q, k, alpha, r0, tt, t_heat, split=False)
        r_t_det = dp_model_recomp_heating_cooling(fs, q, k, tt, t_heat, dT_dual, r0, get_gamma4=False)
    else:
        # only run the algorithm for heating
        dT_dual = compute_delta_T_dual_infinite_line_heating(q, k, alpha, r, tt)
        dT_dual_fixed = compute_delta_T_dual_infinite_line_heating(q, k, alpha, r0, tt)
        r_t_det = dp_model_recomp_heating(fs, q, k, tt, dT_dual, r0, get_gamma4=False)
    r_t_det_diff = r - r_t_det
    r_mm = r / 1.0e-3
    r_det_mm = r_t_det / 1.0e-3
    return dT_dual, dT_dual_fixed, r_det_mm, r_mm, r_t_det_diff, tt


def run_dp_model_inv_test(show_plot=False, run_full=False):
    """
    Run the SP model inverse test to show an example of the reconstruction
    :param  show_plot:          True to show the plot
    :param  run_full:           True to run the reconstruction over the heating and cooling curve
                                False to run the reconstruction over only the heating curve
    :return:
    """

    ################################################
    # MAIN
    ################################################
    q = 45
    k = 5.2
    theta_m = 0.59
    theta_o = 9.2e-3
    theta_w = 0.40
    C = get_volumetric_heat_capacity(theta_m, theta_o, theta_w)
    alpha = get_alpha_k_C(k, C)
    r0 = 6e-3
    fs = 12
    t0 = 0
    max_delta_r = 5.0e-3
    timestep = 1 / fs

    t_heat = 8      # 8 seconds for heating (when working with the full curve)
    if run_full:
        t1 = 60*3   # 3 minutes for heating and cooling
    else:
        t1 = 8      # 8 seconds for heating

    # create the time vector
    t = np.arange(t0, t1, timestep)
    # required since t0 = 0 and this cannot be used since Ei(x=0) is not defined
    rstep = length(t)-1

    ###########################################################
    # LINEAR CHANGE (INCREASE)
    ###########################################################
    print('Computing for linear increase in r(t)...')
    r = r0 + np.linspace(0, max_delta_r, rstep)
    dT_dual, dT_dual_fixed, r_det_mm, r_mm, r_t_det_diff, tt = \
        run_plot_output(alpha, fs, k, q, r0, run_full, t_heat, t, r)

    ###########################################################
    #  LINEAR CHANGE (DECREASE)
    ###########################################################
    print('Computing for linear decrease in r(t)...')
    rdec = r[::-1]
    r0_dec = rdec[0]
    dT_dual_dec, dT_dual_fixed_dec, r_det_mm_dec, r_mm_dec, r_t_det_diff_dec, tt = \
        run_plot_output(alpha, fs, k, q, r0_dec, run_full, t_heat, t, rdec)

    ###########################################################
    # RANDOM WALK
    ###########################################################
    print('Computing for Brownian random walk in r(t)...')
    nr = length(r)
    rb = brownian_noise_norm(nr, r0, r0 + max_delta_r)
    rb0 = rb[0]
    dT_dual_bn, dT_dual_fixed_bn, r_det_mm_bn, r_mm_bn, r_t_det_diff_bn, tt = \
        run_plot_output(alpha, fs, k, q, rb0, run_full, t_heat, t, rb)

    ###########################################################
    # FIGURES
    ###########################################################
    print('DONE COMPUTING, RUNNING FIGURES')
    fig = plt.figure(num=None, figsize=LARGE_SQUARE_FIGSIZE)

    # INCREASE
    # (a)
    ax = fig.add_subplot(3, 3, 1)
    ax.plot(tt, dT_dual_fixed, label='DP Fixed Radius', color='orange')
    ax.plot(tt, dT_dual, label='DP Variable Radius', color="cornflowerblue")
    # ax.set_xlabel('Time (s)')
    ticklabels_off_x()
    ax.set_ylabel(create_label('$\Delta T \hspace{1}$', 'K'))
    ax.legend()
    ax.set_title('(a)', loc='center')

    # (b)
    ax = fig.add_subplot(3, 3, 2)
    ax.plot(tt, r_mm, color=MODELLED_COLOR, linewidth=4, label='Forward')
    ax.plot(tt, r_det_mm, color='grey', ls='dashed', label='Inverse')
    ax.set_ylabel(create_label('$r\hspace{0.3}(\hspace{0.3}t\hspace{0.3})$', 'mm'))
    # ax.set_xlabel('Time (s)')
    ticklabels_off_x()
    ax.legend()
    ax.set_title('(b)', loc='center')

    # (c)
    ax = fig.add_subplot(3, 3, 3)
    ax.plot(tt, r_t_det_diff, color='grey', label='difference')
    ax.set_ylabel(create_label('$r\hspace{0.3}(\hspace{0.3}t\hspace{0.3})\hspace{1}$ Forward - Inverse', 'mm'))
    # ax.set_xlabel('Time (s)')
    ticklabels_off_x()
    ax.set_title('(c)', loc='center')

    # DECREASE
    # (d)
    ax = fig.add_subplot(3, 3, 4)
    ax.plot(tt, dT_dual_fixed_dec, label='DP Fixed Radius', color='orange')
    ax.plot(tt, dT_dual_dec, label='DP Variable Radius', color="cornflowerblue")
    ticklabels_off_x()
    # ax.set_xlabel('Time (s)')
    ax.set_ylabel(create_label('$\Delta T \hspace{1}$', 'K'))
    # ax.legend()
    ax.set_title('(d)', loc='center')

    # (e)
    ax = fig.add_subplot(3, 3, 5)
    ax.plot(tt, r_mm_dec, color=MODELLED_COLOR, linewidth=4, label='Forward')
    ax.plot(tt, r_det_mm_dec, color='grey', ls='dashed', label='Inverse')
    ax.set_ylabel(create_label('$r\hspace{0.3}(\hspace{0.3}t\hspace{0.3})$', 'mm'))
    ticklabels_off_x()
    # ax.set_xlabel('Time (s)')
    # ax.legend()
    ax.set_title('(e)', loc='center')

    # (f)
    ax = fig.add_subplot(3, 3, 6)
    ax.plot(tt, r_t_det_diff_dec, color='grey', label='difference')
    ax.set_ylabel(create_label('$r\hspace{0.3}(\hspace{0.3}t\hspace{0.3})\hspace{1}$ Forward - Inverse', 'mm'))
    ticklabels_off_x()
    # ax.set_xlabel('Time (s)')
    ax.set_title('(f)', loc='center')

    # RANDOM WALK
    # (g)
    ax = fig.add_subplot(3, 3, 7)
    ax.plot(tt, dT_dual_fixed_bn, label='DP Fixed Radius', color='orange')
    ax.plot(tt, dT_dual_bn, label='DP Variable Radius', color="cornflowerblue")
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(create_label('$\Delta T \hspace{1}$', 'K'))
    # ax.legend()
    ax.set_title('(g)', loc='center')

    # (h)
    ax = fig.add_subplot(3, 3, 8)
    ax.plot(tt, r_mm_bn, color=MODELLED_COLOR, linewidth=4, label='Forward')
    ax.plot(tt, r_det_mm_bn, color='grey', ls='dashed', label='Inverse')
    ax.set_ylabel(create_label('$r\hspace{0.3}(\hspace{0.3}t\hspace{0.3})$', 'mm'))
    ax.set_xlabel('Time (s)')
    ax.legend()
    ax.set_title('(h)', loc='center')

    # (i)
    ax = fig.add_subplot(3, 3, 9)
    ax.plot(tt, r_t_det_diff_bn, color='grey', label='difference')
    ax.set_ylabel(create_label('$r\hspace{0.3}(\hspace{0.3}t\hspace{0.3})\hspace{1}$ Forward - Inverse', 'mm'))
    ax.set_xlabel('Time (s)')
    ax.set_title('(i)', loc='center')

    # SHOW AND SAVE FIGURE
    plt.tight_layout()
    plt.savefig(FIG_PATH + DP_SYNTH_FILENAME_EXAMPLE + PLOT_EXTENSION)
    if show_plot:
        plt.show()


def main():
    show_plot = False
    run_dp_model_inv_test(show_plot=show_plot, run_full=True)


if __name__ == '__main__':
    main()

