import numpy as np
from sp_model_signal import sp_model_nominal
from sp_model_signal import sp_model_signal_inv, dT_sp_model_signal
from gen_time_vec import gen_time_vec
from cut_to_zero import cut_to_zero
from constant_labels import create_label
from inverse_model_nominal import obtain_sp_vars_from_curve
from ticklabels_off import ticklabels_off_x
from constants import *

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('ggplot')
mpl.rcParams["mathtext.fontset"] = "stix"


def run_sp_model_example(q, k, t, b, c, d, show_plot=False):
    """
    Run the SP model example and the inverse to show how a synthetic curve can be obtained
    :param q:
    :param k:
    :param t:
    :param b:
    :param c:
    :param d:
    :return:
    """
    dT = sp_model_nominal(q, k, t, b, c, d)
    dT_synth0, t0 = cut_to_zero(dT, t)
    dt = t0[1] - t0[0]
    kdet, bdet, tth, dT_sig_pro = sp_model_signal_inv(dT_synth0, t0, dt, q, output_model_array=True,
                                                      cutoff=0.5, filter=False)
    dT_sig_pro_comp = dT_sp_model_signal(q, kdet, tth, bdet)
    dT_sig_pro_comp_diff = dT_sig_pro - dT_sig_pro_comp
    kd, bd, cd, dd = obtain_sp_vars_from_curve(q, t0, dT_synth0, kdet)
    dT_comp = sp_model_nominal(q, kd, t0, bd, cd, dd)
    dT_comp_diff = dT_synth0 - dT_comp

    fig = plt.figure(num=None, figsize=SQUARE_FIGSIZE)

    # (a)
    ax = fig.add_subplot(2, 2, 1)
    # ax.set_xlabel('Time (s)')
    ticklabels_off_x()
    ax.set_ylabel(create_label('$\Delta T$', 'K'))
    ax.plot(t0, dT_synth0, color=MODELLED_COLOR, linewidth=4, label='Forward')
    ax.plot(t0, dT_comp, color='grey', ls='dashed', label='Inverse')
    ax.legend()
    ax.set_title('(a)', loc='center')

    # (b)
    ax = fig.add_subplot(2, 2, 2)
    # ax.set_xlabel('Time (s)')
    ticklabels_off_x()
    ax.set_ylabel(create_label('$\Delta T \hspace{1}$ Forward - Inverse', 'K'))
    ax.plot(t0, dT_comp_diff, color=GREY_COLOR, label='Difference', ls='solid')
    ax.legend()
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax.set_title('(b)', loc='center')

    # (c)
    ax = fig.add_subplot(2, 2, 3)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(create_label('$\Delta \Gamma _{5} \hspace{0.3} (\hspace{0.3} t \hspace{0.3})$', ''))
    ax.plot(tth, dT_sig_pro_comp, color=MODELLED_COLOR, linewidth=4, label='Forward (Processed)')
    ax.plot(tth, dT_sig_pro, color='grey', ls='dashed', label='Theoretical')
    ax.set_title('(c)', loc='center')
    ax.legend()

    # (d)
    ax = fig.add_subplot(2, 2, 4)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(create_label('$\Delta \Gamma _{5} \hspace{0.3} (\hspace{0.3} t \hspace{0.3}) \hspace{1}$ '
                               'Forward - Inverse', ''))
    ax.plot(tth, dT_sig_pro_comp_diff, color=GREY_COLOR, ls='solid')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax.set_title('(d)', loc='center')

    plt.tight_layout()
    plt.savefig(FIG_PATH + SP_SYNTH_FILENAME_EXAMPLE + PLOT_EXTENSION)
    if show_plot:
        plt.show()


def run_sp_synth(show_plot=False):
    T = 8.0
    fs = FS_SAMPLE
    dt = 1 / fs
    t = gen_time_vec(fs, T)
    t = t[1:]   # time vector cannot start at zero

    # example values
    q = 45
    k = 5.2
    b = 0.01642178907454129
    c = 0.20330724016048247
    d = 0.4026793334317176
    run_sp_model_example(q, k, t, b, c, d, show_plot)


def main():
    run_sp_synth(show_plot=False)


if __name__ == '__main__':
    main()
