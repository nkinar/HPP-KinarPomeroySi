from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
from latex_float import latex_float
import numpy as np

"""
NOTES on offset API
1. https://matplotlib.org/api/offsetbox_api.html
"""


def place_rmse_mb_pd_on_plot(rmse, mb, pd, units, dec_places, loc, is_scientific, use_pd=False):
    """
    Function to place the RMSE and MB on the most recent plot
    :param rmse:            the rmse as a number
    :param mb:              as the mean bias as a number
    :param units:           units for mb
    :param dec_places:      number of decimal places for the mean bias
    :param loc:             default location to place the strings
    :param is_scientific:   True to use scientific notation
    :return:
    """
    ff = '{:.' + str(dec_places)
    if is_scientific:
        ff += 'e}'
    else:
        ff += '}'
    frmse = ff.format(rmse)
    fmb = ff.format(mb)
    frmse_format = latex_float(frmse)
    fmb_format = latex_float(fmb)
    s = r'$\mathrm{RMSD} = ' + frmse_format + r'\hspace{0.3} \mathrm{' + units + r'}$' + '\n' + \
        r'$\mathrm{MB} = ' + fmb_format + r'\hspace{0.3} \mathrm{' + units + r'}$'
    if use_pd:
        if pd > 1:
            fpd = str(int(np.round(pd)))
        else:
            fpd = latex_float(ff.format(pd))
        s += '\n'
        s +=  r'$\mathrm{PD} = ' + fpd + r'\hspace{0.3} \mathrm{' + '\%' + r'}$'
    ax = plt.gca()
    anchored_text = AnchoredText(s, loc)
    ax.add_artist(anchored_text)
    # DONE



