import numpy as np

def get_max(y):
    """
    Helper function to obtain the index and the value from an array
    :param y:
    :return:
    """
    idx = np.argmax(y)
    val = y[idx]
    return idx, val


def get_min(y):
    """
    Helper function to obtain the index and the value from an array
    :param y:
    :return:
    """
    idx = np.argmin(y)
    val = y[idx]
    return idx, val
