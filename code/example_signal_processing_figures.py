from recomp_theta_rho import recomp_theta_rho_additional
from constant_labels import create_label
from constants import *
from set_xlim_linear import set_xlim_linear
from prettytable import PrettyTable
from comparisons import compute_percentage_diff, compute_mb, compute_rmse
from float_round import float_round

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('ggplot')
mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams.update({'font.size': FONT_SIZE})


def run_example_signal_processing_figure(theta_known, density_known, path, start_string, additional_text,
                                         q, t_heat, t_total, end_trim_time, fname,
                                         show_plot=False):
    """
    Run the example signal processing to construct a figure
    :param path:
    :param start_string:
    :param additional_text:
    :param q:
    :param t_heat:
    :param t_total:
    :param end_trim_time:
    :param fname:
    :return:
    """
    theta_w_nom, rho_nom, theta_w_sig, rho_sig, \
    t2_trim1, delta_T2_trim1, dT_synth_nom, dT_synth_sig, \
    cut_time, t2_trim1_cut, r_t_heating_cooling, \
    alpha_det, alpha_sig, k_det, kdet_sig, \
    t1_trim_heating1, delta_T1_trim_heating1, dT_synth_sp, \
    delta_T1_trim, t1_trim, \
    tlinear, kdet_linear, qav = \
        recomp_theta_rho_additional(path, start_string, additional_text, q, t_heat, t_total, end_trim_time)

    fig = plt.figure(num=None, figsize=HLONG_FIGSIZE_M)

    # (a) Single Probe
    ax = fig.add_subplot(1, 3, 1)
    ax.plot(t1_trim, delta_T1_trim, label='Measured', color=MEASURED_COLOR)
    ax.plot(t1_trim_heating1, dT_synth_sp, label='Modelled', color=MODELLED_COLOR)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(create_label('$\Delta T \hspace{1}$', 'K'))
    ylim = list(ax.get_ylim())
    ylim[1] += 0.20*(ylim[1]-ylim[0])
    ax.set_ylim(ylim)
    ax.legend(loc='best')
    ax.set_title('(a)', loc='left')

    # (b) Dual Probe
    ax = fig.add_subplot(1, 3, 2)
    ax.plot(t2_trim1, delta_T2_trim1, label='Measured', color=MEASURED_COLOR)
    ax.plot(t2_trim1, dT_synth_nom, label='Modelled', color=MODELLED_COLOR)
    ax.axvline(cut_time, linestyle='--', label='Peak Time', color=MODELLED_COLOR)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(create_label('$\Delta T \hspace{1}$', 'K'))
    ylim = list(ax.get_ylim())
    ylim[1] += 0.20*(ylim[1]-ylim[0])
    ax.set_ylim(ylim)
    ax.legend(loc='best')
    ax.set_title('(b)', loc='left')

    # (d) r(t)
    ax = fig.add_subplot(1, 3, 3)
    ax.plot(t2_trim1_cut, r_t_heating_cooling / 1.0e-3, label='Effective Radius', color='grey')
    ax.legend(loc='upper left')
    ax.set_ylabel(create_label('$r\hspace{0.3}(\hspace{0.3}t\hspace{0.3})$', 'mm'))
    ax.set_xlabel('Time (s)')
    ylim = list(ax.get_ylim())
    ylim[1] += 0.20*(ylim[1]-ylim[0])
    ax.set_ylim(ylim)
    ax.set_title('(c)', loc='left')

    plt.tight_layout()
    plt.savefig(FIG_PATH + fname + PLOT_EXTENSION)

    # compute percentage differences
    pd_theta_sig = float_round(compute_percentage_diff(theta_known, theta_w_sig), NOM_DP)
    pd_theta_nom = float_round(compute_percentage_diff(theta_known, theta_w_nom), NOM_DP)
    diff_theta_sig = float_round(theta_known - theta_w_sig, NOM_DP)
    diff_theta_nom = float_round(theta_known - theta_w_nom, NOM_DP)

    pd_rho_sig = float_round(compute_percentage_diff(density_known, rho_sig), NOM_DP)
    pd_rho_nom = float_round(compute_percentage_diff(density_known, rho_nom), NOM_DP)
    diff_rho_sig = float_round(density_known - rho_sig, NOM_DP)
    diff_rho_nom = float_round(density_known - rho_nom, NOM_DP)

    # create HTML table to show values
    tab = PrettyTable()
    tab.field_names = ["Type", "Signal Hybrid SP and DP", "DP nominal", "late-time SP nominal"]
    tab.add_row(["k", float_round(kdet_sig, NOM_DP), float_round(k_det, NOM_DP), float_round(kdet_linear, NOM_DP)])
    tab.add_row(["alpha", float_round(alpha_sig, NOM_DP), float_round(alpha_det, NOM_DP), None])
    tab.add_row(["theta_w", float_round(theta_w_sig, NOM_DP), float_round(theta_w_nom, NOM_DP), None])
    tab.add_row(["rho", float_round(rho_sig, NOM_DP), float_round(rho_nom, NOM_DP), None])
    print(tab)

    # create HTML table with differences
    tabd = PrettyTable()
    tabd.field_names = ["Type", "PD Signal (%)", "Diff signal", "PD nom (%)", "Diff nom"]
    tabd.add_row(["theta_w", float_round(pd_theta_sig, NOM_DP),
                  float_round(diff_theta_sig, NOM_DP),
                  float_round(pd_theta_nom, NOM_DP),
                  float_round(diff_theta_nom, NOM_DP)])
    tabd.add_row(["rho",  float_round(pd_rho_sig, NOM_DP),
                  float_round(diff_rho_sig, NOM_DP),
                  float_round(pd_rho_nom, NOM_DP),
                  float_round(diff_rho_nom, NOM_DP)])
    print(tabd)

    # save out the tables to HTML file
    print('Saving tables....')
    with open(TABLE_PATH + fname + '-table1' + HTML_EXT, 'w') as f:
        f.write(tab.get_html_string())
        f.close()

    with open(TABLE_PATH + fname + '-table2' + HTML_EXT, 'w') as f:
        f.write(tabd.get_html_string())
        f.close()
    print('DONE saving tables...')

    if show_plot:
        plt.show()


def run_example_signal_sand(show_plot=False):
    path = SAND_PATH_EXAMPLE
    start_string = 'sand-'
    additional_text = '-rep3'
    q = 45
    t_heat = 8
    t_total = 3 * SECONDS_IN_MIN
    end_trim_time = 0
    fname = SAND_EXAMPLE
    run_example_signal_processing_figure(THETA_NOMINAL_SAND, DENSITY_NOMINAL_SAND, path,
                                         start_string, additional_text, q, t_heat, t_total, end_trim_time, fname,
                                         show_plot=show_plot)


def run_example_signal_peat(show_plot=False):
    path = PEAT_PATH_EXAMPLE
    start_string = 'peat-'
    additional_text = '-rep4'
    q = 20
    t_heat = 89
    t_total = 3 * SECONDS_IN_MIN
    end_trim_time = 0
    fname = PEAT_EXAMPLE
    run_example_signal_processing_figure(THETA_NOMINAL_PEAT, DENSITY_NOMINAL_PEAT, path,
                                         start_string, additional_text, q, t_heat, t_total, end_trim_time, fname,
                                         show_plot=show_plot)


def main():
    print('Calculating for sand example...')
    run_example_signal_sand(show_plot=False)
    print('Calculating for peat example...')
    run_example_signal_peat(show_plot=False)


if __name__ == '__main__':
    main()

