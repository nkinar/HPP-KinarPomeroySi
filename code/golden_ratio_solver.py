import numpy as np

def golden_ratio_solver(f, ain, bin, max_iter=100000, tol=1.0e-15):
    """
    Golden ratio solver to minimize f as a function
    :param f:           as the function to minimize f(x) = 0
    :param ain:         as the left starting point of the interval
    :param bin:         as the right starting point of the interval
    :return:
    """
    gr = (np.sqrt(5.0) + 1.0) * 0.5
    a = ain
    b = bin
    for k in range(max_iter):
        c = b - (b - a) / gr
        d = a + (b - a) / gr
        if np.abs(c-d) < tol:
            break
        if f(c) < f(d):
            b = d
        else:
            a = c
    out = 0.5*(b + a)
    return out


############################################################################
# EXAMPLE CODE
############################################################################


def main():
    def _f(x):
        out = (x-1)**2
        return out
    a = -2
    b = 3
    out = golden_ratio_solver(_f, a, b)
    print('out = ', out) # min should be at 1.0


if __name__ == '__main__':
    main()

