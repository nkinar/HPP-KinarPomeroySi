import matplotlib.pyplot as plt

def set_x_scientific():
    ax = plt.gcf()
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

def set_y_scientific():
    ax = plt.gca()
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))