import numpy as np
from get_size import length

"""
This is a collection of routines used to compute the forward and inverse
2-point derivative. 

Tikhonov regularization is also used here as well.
"""


def inverse_first_derivative(xprime, dt, a0):
    """
    Take the inverse of the first derivative using recursion
    :param xprime:          as the derivative vector
    :param dt:              as the timestep between vectors
    :param a0:              as the known first element (boundary condition)
    :return:
    """
    n = length(xprime)
    n1 = n + 1
    out = np.zeros(n1)
    out[0] = a0
    cnt = 0
    for i in range(1, n1):
        out[i] = xprime[cnt] * dt + out[cnt]
        cnt += 1
    return out


def forward_first_derivative(x, dx):
    """
    Take the first derivative using a two-point stencil
    :param x:       as the vector
    :param dx:      as the step between elements
    :return:
    """
    n = length(x)
    nsub1 = n - 1
    out = np.zeros(nsub1)
    cnt = 0
    for i in range(1, n):
        out[cnt] = (x[i] - x[i-1])/dx
        cnt += 1
    return out

