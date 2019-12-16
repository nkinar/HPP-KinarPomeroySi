import numpy as np
from scipy.signal import butter, lfilter, freqz, filtfilt


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5, both=True):
    b, a = butter_lowpass(cutoff, fs, order=order)
    if both:
        y = filtfilt(b, a, data)
    else:
        y = lfilter(b, a, data)
    return y