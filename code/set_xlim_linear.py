import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def set_xlim_linear(start, end, num, integer=True):
    """
    Set the xlim as a linear axis
    :param start:   first element
    :param end:     last element
    :param num:     number of ticks
    :return:
    """
    plt.xlim([start, end])
    tt = np.linspace(start, end, num)
    ax = plt.gca()
    ax.set_xticks(tt)
    if integer:
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))