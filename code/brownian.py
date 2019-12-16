import numpy as np
from normalize import normalize

def white_noise(n):
    """
    Generate white noise as numbers from a normal distribution
    The power spectrum is proportional to P(f) = 1/f^beta = 1
    where beta = 0
    :param n:               as the length of the white noise
    :return:                white noise vector
    """
    mean = 0
    std = 1
    out = np.random.normal(mean, std, size=n)
    return out


def brownian_noise(n):
    """
    Generate Brownian noise by integration of white noise
    The power spectrum is proportional to P(f) = 1/f^beta = 1/f^2
    where beta = 2
    :param n:           as the length of the Brownian noise
    :return:            Brownian noise vector
    """
    white_vec = white_noise(n)
    out = np.cumsum(white_vec)
    return out


def brownian_noise_norm(n, a, b):
    """
    Generate Brownian noise that is normalized to an interval
    :param n:   as the number of points
    :param a:   as the left endpoint of the interval (min)
    :param b:   as the right endpoint of the interval (max)
    :return:
    """
    bn = brownian_noise(n)
    noise = normalize(bn, a, b)
    return noise

####################################################################################

def main():
    """
    Sample test code
    :return:
    """
    import matplotlib.pyplot as plt
    block = True
    n = 1000
    a = 5.0e-3
    b = 12e-3
    noise = brownian_noise_norm(n, a, b)
    plt.figure()
    plt.plot(noise)
    plt.show(block)


if __name__ == '__main__':
    main()


