import numpy as np
import pickle
from constant_labels import create_label
from get_size import length
from set_axis_integer import set_x_axis_integer
from ticklabels_off import ticklabels_off_y, ticklabels_off_x
from format_scientific_axis import set_y_scientific
from constants import *

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('ggplot')
mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams.update({'font.size': FONT_SIZE})


def plot_main_all_data(show_plot=False):
    print('Loading data...')
    experiments = pickle.load(open(DATA_PICKLE_FILE_VEC, "rb"))

    theta_nom_sand, rho_nom_sand, theta_sig_sand, rho_sig_sand, \
        alpha_det_sand, alpha_sig_sand, k_det_sand, kdet_sig_sand, kdet_linear_sand, \
    qav_vec_sand, qknown_vec_sand = experiments[SAND_NAME]

    theta_nom_peat, rho_nom_peat, theta_sig_peat, rho_sig_peat, \
        alpha_det_peat, alpha_sig_peat, k_det_peat, kdet_sig_peat, kdet_linear_peat, \
    qav_vec_peat, qknown_vec_peat \
        = experiments[PEAT_NAME]
    print('Done loading data...')

    # first sequence of experiments
    ns = length(theta_nom_sand)
    xs = np.linspace(1, ns, ns)

    # second sequence of experiments
    nnp = length(theta_nom_peat)
    xp = np.linspace(xs[-1]+1, xs[-1]+nnp, nnp)

    #################################################################

    fig = plt.figure(num=None, figsize=SQUARE_FIGSIZE)

    lim_theta = [0, 1.0]
    lim_rho = [400, 2500]

    # (a) theta_sand
    ax = fig.add_subplot(2, 2, 1)
    ax.plot(xs, theta_nom_sand, 's', label='Heating and Cooling DP')
    ax.plot(xs, theta_sig_sand, 'o', label='Signal Processing SP and DP')
    ax.axhline(THETA_NOMINAL_SAND, linestyle='--', label='Nominal Value', color='gray')
    ax.set_ylim(lim_theta)
    set_x_axis_integer()
    ax.set_ylabel(create_label(r'$\theta_w$', ''))
    ticklabels_off_x()
    ax.legend()
    ax.set_title('(a)', loc='left')

    # (b) theta_peat
    ax = fig.add_subplot(2, 2, 2)
    ax.plot(xp, theta_nom_peat, 's', label='Heating and Cooling DP')
    ax.plot(xp, theta_sig_peat, 'o', label='Signal Processing SP and DP')
    ax.axhline(THETA_NOMINAL_PEAT, linestyle='--', label='Nominal Value', color='gray')
    ax.set_ylim(lim_theta)
    ticklabels_off_y()
    ticklabels_off_x()
    set_x_axis_integer()
    ax.set_title('(b)', loc='left')

    # (c) density of sand
    ax = fig.add_subplot(2, 2, 3)
    ax.plot(xs, rho_nom_sand, 's', label='Heating and Cooling DP')
    ax.plot(xs, rho_sig_sand, 'o', label='Signal Processing SP and DP')
    ax.axhline(DENSITY_NOMINAL_SAND, linestyle='--', label='Nominal Value', color='gray')
    ax.set_ylim(lim_rho)
    set_x_axis_integer()
    ax.set_ylabel(create_label(r'$\rho$', 'kg m^-3'))
    ax.set_xlabel('Experiment #')
    ax.set_title('(c)', loc='left')

    # (d) density of peat
    ax = fig.add_subplot(2, 2, 4)
    ax.plot(xp, rho_nom_peat, 's', label='Heating and Cooling DP')
    ax.plot(xp, rho_sig_peat, 'o', label='Signal Processing SP and DP')
    ax.axhline(DENSITY_NOMINAL_PEAT, linestyle='--', label='Nominal Value', color='gray')
    ax.set_ylim(lim_rho)
    ticklabels_off_y()
    set_x_axis_integer()
    ax.set_xlabel('Experiment #')
    ax.set_title('(d)', loc='left')

    plt.tight_layout()
    plt.savefig(FIG_PATH + PLOTS_ALL_THETA_RHO_FILENAME + PLOT_EXTENSION)

    if show_plot:
        block = False
        plt.show(block)

    ###########################################################################

    lim_k = [0, 9]
    lim_a = [1e-7, 6e-6]

    fig = plt.figure(num=None, figsize=SQUARE_FIGSIZE)

    # (a) k for sand
    ax = fig.add_subplot(2, 2, 1)
    ax.plot(xs, k_det_sand, 's', label='Heating and Cooling DP')
    ax.plot(xs, kdet_sig_sand, 'o', label='Signal Processing SP')
    ax.plot(xs, kdet_linear_sand, '^', label='Late-Time SP')
    ax.set_ylim(lim_k)
    set_x_axis_integer()
    ax.set_ylabel(create_label(r'$k$', 'W m^-1 K^-1'))
    ticklabels_off_x()
    ax.legend()
    ax.set_title('(a)', loc='center')

    # (b) alpha for sand
    ax = fig.add_subplot(2, 2, 3)
    ax.plot(xs, alpha_det_sand, 's', label='Heating and Cooling DP')
    ax.plot(xs, alpha_sig_sand, 'o', label='Signal Processing SP and DP')
    ax.set_ylim(lim_a)
    set_x_axis_integer()
    set_y_scientific()
    ax.set_ylabel(create_label(r'$\alpha$', 'm^2 s^-1'))
    ax.set_xlabel('Experiment #')
    ax.set_title('(c)', loc='center')
    ax.legend()

    # (c) k for peat
    ax = fig.add_subplot(2, 2, 2)
    ax.plot(xp, k_det_peat, 's', label='Heating and Cooling DP')
    ax.plot(xp, kdet_sig_peat, 'o', label='Signal Processing SP')
    ax.plot(xp, kdet_linear_peat, '^', label='Late-Time SP')
    ax.set_ylim(lim_k)
    set_x_axis_integer()
    ticklabels_off_x()
    ticklabels_off_y()
    ax.set_title('(b)', loc='center')

    # (d) alpha for peat
    ax = fig.add_subplot(2, 2, 4)
    ax.plot(xp, alpha_det_peat, 's', label='Heating and Cooling DP')
    ax.plot(xp, alpha_sig_peat, 'o', label='Signal Processing SP')
    set_x_axis_integer()
    set_y_scientific()
    ax.set_ylim(lim_a)
    ticklabels_off_y()
    ax.set_xlabel('Experiment #')
    ax.set_title('(d)', loc='center')

    plt.tight_layout()
    plt.savefig(FIG_PATH + PLOTS_ALL_THERMAL_FILENAME + PLOT_EXTENSION)

    if show_plot:
        block = True
        plt.show(block)
    # DONE


def main():
    plot_main_all_data(show_plot=False)


if __name__ == '__main__':
    main()

