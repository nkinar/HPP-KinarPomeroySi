import numpy as np

def running_mean(x, N):
    out = np.convolve(x, np.ones((N,))/N, mode='same')
    return out