import numpy as np
from get_size import length


def _first(x, k, dx2):
    out = (x[k+1] - x[k-1])/dx2
    return out


def _second(x, k, dx2):
    out = (x[k - 1] - 2 * x[k] + x[k + 1]) / dx2
    return out


def _deriv_3pt_common(x, dx, t, typ):
    if typ == 'first':
        f = _first
        dx2 = 2.0 * dx
    elif typ == 'second':
        f = _second
        dx2 = dx**2
    else:
        raise ValueError('Type of derivative not recognized')
    n = length(x)
    nn = n - 2
    out = np.zeros(nn)
    cnt = 0
    tt = t[1:-1]
    for k in range(1, n-1):
        out[cnt] = f(x, k, dx2)
        cnt += 1
    return tt, out


def first_derivative_3pt(x, dx, t):
    """
    Compute the first derivative using 3 pt stencil
    :param x:       as the array to compute the derivative
    :param dx:      as the time or spatial step
    :param t:       as the time vector to be returned as a trimmed vector
    :return:        (tt, out)
    tt              as the truncated time vector
    out             as the first derivative computed using a 3-point stencil
    """
    return  _deriv_3pt_common(x, dx, t, 'first')


def second_derivative_3pt(x, dx, t):
    """
    Compute the second derivative using three-point stencil
    :param x:
    :param dx:
    :param t:
    :return:
    """
    return  _deriv_3pt_common(x, dx, t, 'second')


#######################################################
# EXAMPLE USAGE
#######################################################

def main():
    import matplotlib.pyplot as plt
    t = np.linspace(0, 10, 100)
    dx = t[1] - t[0]
    x = 3*t + t**2

    tt, dy = first_derivative_3pt(x, dx, t)
    dy_comp = 3 + 2*tt
    diff = dy-dy_comp

    block = False
    plt.figure()
    plt.plot(t, x)
    plt.show(block)

    block = False
    plt.figure()
    plt.plot(tt, dy)
    plt.plot(tt, dy_comp)
    plt.show(block)

    block = True
    plt.figure()
    plt.plot(tt, diff)
    plt.show(block)



if __name__ == '__main__':
    main()


