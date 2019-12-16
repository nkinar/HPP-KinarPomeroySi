from get_size import length
import numpy as np

def safe_div(x, y, replace=None):
    if replace is None:
        replace = float('NaN')
    n = length(x)
    out = np.zeros(n)
    if n != length(y):
        raise ValueError('save_div: cannot divide')
    for k in range(n):
        xx = x[k]
        yy = y[k]
        if yy == 0:
            out[k] = replace
        else:
            out[k] = xx/yy
    return out