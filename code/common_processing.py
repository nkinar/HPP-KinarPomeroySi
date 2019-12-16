import numpy as np
from constants import *


def get_total_time_signal(n, fs):
    """
    Get the total time of the signal
    :param n:       as the length of the signal
    :param fs:      as the sampling freqency
    :return:
    """
    t_total = float(n)*fs
    return t_total


def get_dt(fs):
    """
    Obtain the sampling step
    :param fs:
    :return:
    """
    return 1/fs


def av_first_part_of_curve(T, fs, tav):
    """
    Average the first part of the curve to obtain the initial temperature
    :param T:           as the temperature curve            [K or C]
    :param fs:          as the sampling rate                [Hz]
    :param tav:         as the time over which to sample    [s]
    :return:
    """
    nav = int(np.ceil(tav*fs))
    out = np.average(T[0:nav])
    return out


def compute_difference_curve(T, fs, tav):
    """
    Obtain a difference curve from T
    :param T:       as the temperature curve [K or C]
    :param fs:      as the sampling rate [Hz]
    :param tav:     as the time over which to sample [s]
    :return:
    """
    av = av_first_part_of_curve(T, fs, tav)
    Tdiff = T-av
    return Tdiff



