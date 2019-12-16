import numpy as np

def gen_time_vec(fs, T, tstart=0):
    """
    Generate a time vector given the sampling rate and total time
    :param fs:      as the sampling frequency (Hz)
    :param T:       as the total time (s)
    :param tstart:  as the starting time (s)  [normally zero]
    :return:
    """
    n = int(np.ceil(T*fs))
    dt = 1.0 / fs
    tout = np.linspace(tstart, T-dt, n)
    return tout