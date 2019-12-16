import matplotlib.pyplot as plt

def ticklabels_off_x():
    """
    Turn the ticklabels off on the y axis
    :return:
    """
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])


def ticklabels_off_y():
    """
    Turn the ticklabels off on the x axis
    :return:
    """
    frame1 = plt.gca()
    frame1.axes.yaxis.set_ticklabels([])