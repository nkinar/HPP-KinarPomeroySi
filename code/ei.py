import numpy as np
from scipy.special import expi
from OptimHelper import OptimHelper
from golden_ratio_solver import golden_ratio_solver
from float_comparison import get_machine_epsilon
from get_size import length


def ei(x):
    """
    Compute the exponential integral
    :param x:   as the argument of the function
    :return:    output from the function
    """
    return expi(x)

def e1(x):
    """
    Compute the exponential integral e1(x)
    :param x:
    :return:
    """
    return -ei(-x)


def ei_inv_neg_array(known):
    """
    Take the inverse of Ei(x) with x < 0
    :param known:
    :return:
    """
    n = length(known)
    out = np.zeros(n)
    for k in range(n):
        out[k] = ei_inv_neg(known[k])
    return out

def ei_inv_neg(known):
    """
    Take the inverse of Ei(x) with x < 0
    Since the function has a discontinuity at x = 0, there are two branches.
    :param known:       as the known function output
    :param is_neg:      True if the original x is negative
    :return:
    """
    eps = get_machine_epsilon()
    oh = OptimHelper(ei, known)
    limit = 1.0/(known**2)
    a = -limit
    b = eps
    out = golden_ratio_solver(oh.fcall, a, b, tol=1.0e-30)
    return out


############################################################################
# EXAMPLE CODE
############################################################################


def main():
    x = -1.2
    first = ei(x)
    print('first = ', first)
    inv_first = ei_inv_neg(first)
    print('inv_first = ', inv_first)



if __name__ == '__main__':
    main()


