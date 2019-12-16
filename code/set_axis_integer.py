from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt

def set_x_axis_integer():
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))


def set_y_axis_integer():
    ax = plt.figure().gca()
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
