import numpy as np


def get_eps():
    """
    Similar to the eps(1) function in Matlab
    """
    eps = np.spacing(1.0)
    return eps