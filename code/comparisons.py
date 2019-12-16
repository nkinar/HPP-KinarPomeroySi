import numpy as np
from get_size import length
from safe_div import safe_div


def compute_rmse(known, model):
    """
    Computes the RMSE (Root Mean Squared Error) and returns this in the same units as the data
    :param known:
    :param model:
    :return:
    """
    k = np.asarray(known)
    m = np.asarray(model)
    n = k.size
    term = np.sqrt(np.sum((m-k)**2)/float(n))
    return term


def compute_mb_same_units(known, model):
    """
    Computes the Mean Bias (MB) and returns this in the same units as the data
    :param known:
    :param model:
    :return:
    """
    k = np.asarray(known)
    m = np.asarray(model)
    n = k.size
    term = np.sum(m-k)/float(n)
    return term


def compute_mb(known, model):
    """
    Computes the Mean Bias (MB)
    :param known:
    :param model:
    :return:
    """
    return compute_mb_same_units(known, model)


def compute_percentage_diff(known, model):
    """
    Compute the percentage difference
    :param known:   as the known value
    :param model:   as the modelled value
    :return:
    """
    k = np.asarray(known)
    m = np.asarray(model)
    term = np.mean((k-m)/k)
    return 100*term


def compute_pd(known, model):
    """
    Convenience function
    :param known:
    :param model:
    :return:
    """
    return compute_percentage_diff(known, model)


def compute_variance(x):
    """
    Compute variance
    :param x:
    :return:
    """
    return np.var(x)


def compute_sd(x):
    """
    Compute the standard deviation
    :param x:
    :return:
    """
    return np.std(x)

