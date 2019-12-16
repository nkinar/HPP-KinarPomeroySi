import numpy as np
from get_size import length


def average_downsample(fsample, fwanted, v):
    """
    Downsample a signal over a series of windows
    :param fsample:     as the actual sampling rate
    :param fwanted:     as the wanted sampling rate
    :param v:
    :return:
    """
    n = int(np.ceil(fsample/fwanted))  # number of points in the window
    nlen = length(v)
    nwin = int(nlen/n)                 # number of windows
    begin = 0
    end = n
    out = np.zeros(nwin)               # the output has the same length as the number of windows
    for k in range(nwin):
        win = v[begin:end]
        av = np.average(win)
        out[k] = av
        begin = end
        end += n
    return out

