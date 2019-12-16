import numpy as np
from gaussian_curve import gaussian
from get_size import length


def gen_radius(r0, r1, t, th, sd, heat_cool=False):
    """
    Generate a radius suitable for demonstrating the inverse
    :param r0:              as the starting radius
    :param r1:              as the ending radius
    :param t:               as the total time vector
    :param th:              as the time of heating
    :param sd:              as the standard deviation
    :param heat_cool:       True to generate heating and cooling curves
    :return:
    """
    n = length(t)               # obtain the length of t
    if heat_cool:               # run over the heating and cooling curve
        a = r1                  # height of the peak at the maximum radius
        b = th                  # height of the peak time
        c = sd                  # standard deviation
        rout = r0 + gaussian(a, b, c, t)
    else:
        rout = np.linspace(r0, r1, n)
    return rout

##############################################################################

def main():
    from gen_time_vec import gen_time_vec
    import matplotlib.pyplot as plt

    heat_cool = True
    fs = 120
    r0 = 6e-3
    r1 = 11e-3
    th = 8
    if heat_cool:
        T = 60
    else:
        T = th
    t = gen_time_vec(fs, T)
    sd = 5
    r = gen_radius(r0, r1, t, th, sd, heat_cool)

    block = True
    plt.figure()
    plt.plot(t, r)
    plt.show(block)


if __name__ == '__main__':
    main()

