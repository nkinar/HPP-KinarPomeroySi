import numpy as np

def gen_mid_vec(rmin, rmid, rmax, n):
    """
    Generate a vector where there is a min and max and a midpoint and there are
    n points distributed between the mid and the closest end
    :param rmin:        as the min point
    :param rmid:        as the mid point
    :param rmax:        as the max point
    :param n:           as the number
    :return:
    """
    initial = np.linspace(rmin, rmid, n)
    last = np.linspace(rmid, rmax, n)[1:]
    full = np.concatenate([initial, last])
    return full

######################################################################################################

def main():
    rmin = 1.0e-3
    rmid = 6.0e-3
    rmax = 11e-3
    n = 10
    out = gen_mid_vec(rmin, rmid, rmax, n)
    print(out)


if __name__ == '__main__':
    main()
