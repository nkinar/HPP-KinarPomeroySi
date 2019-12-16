import numpy as np
from tinv import tinv_two_tail, tinv_one_tail
from get_size import length
from even_odd import is_even
from float_comparison import close
from running_mean_filter import running_mean
from zero_runs import zero_runs

# Indicates what is true or false for changes in the curve
CURVE_CHANGE_TRUE_VAL = 1
CURVE_CHANGE_FALSE_VAL = 0


def obtain_t(x, y):
    """
    Compute the Student's t value between two datasets
    :param x:
    :param y:
    :return: (t, df)
    t = Student's t value
    df = degrees of freedom
    """
    nx = length(x)
    ny = length(y)
    if nx == 0 or ny == 0:
        raise ValueError('obtain_t: datasets cannot have zero elements')
    xbar = np.mean(x)
    ybar = np.mean(y)
    numerator = np.abs(xbar-ybar)
    sx = ((np.sum(x**2)/nx) - (xbar**2))/(nx-1.0)
    sy = ((np.sum(y**2)/ny) - (ybar**2))/(ny-1.0)
    if sx < 0:
        sx = 0
    if sy < 0:
        sy = 0
    denominator = np.sqrt(sx + sy)
    if denominator == 0.0:
        t = 0
    else:
        t = numerator/denominator
    df_numerator = (sx + sy)**2
    df_denominator = ((sx**2)/(nx-1.0)) + ((sy**2)/(ny-1.0))
    if df_denominator == 0.0:
        df = 0.0
    else:
        df = numerator/denominator
    return t, df


class StudentTest:
    def __init__(self, p, two_tail=False):
        if p < 0 or p > 1:
            raise ValueError('StudentTest: p must be greater than 0 and less than 1')
        self.p = p
        self.tinv = tinv_one_tail
        if two_tail:
            self.tinv = tinv_two_tail
        self.t_last = None
        self.df_last = None

    def run(self, x, y):
        t, df = obtain_t(x, y)
        if t == 0 or df == 0:
            return False
        compare = self.tinv(self.p, df)
        if t < compare:
            return False    # null hypothesis cannot be rejected
        return True         # reject the null hypothesis
# DONE


class PeakDetectStudent:
    def __init__(self, p, ws, two_tail=True):
        self.st = StudentTest(p, two_tail)
        if is_even(ws):
            raise ValueError('PeakDetectStudent: Window size must be odd')
        self.ws = ws


    def find_curve_changes(self, x):
        nx = length(x)
        n = nx - 2*self.ws
        start = 0
        end = self.ws
        mid = int(np.floor(0.5 * (start + end)))
        out = np.zeros(nx)  # output is same size as the initial input
        for k in range(n):
            first = x[start:mid]
            second = x[mid+1:end]
            mid += 1
            if first.size == 0 or second.size == 0:
                break
            test = self.st.run(first, second)
            if test:
                out[mid] = CURVE_CHANGE_TRUE_VAL
            else:
                out[mid] = CURVE_CHANGE_FALSE_VAL
            # switch out the endpoints
            start += 1
            end += 1
        return out
# DONE

class FindPlateau:
    def __init__(self, p, ws, two_tail=True):
        self.ws = ws
        self.pds = PeakDetectStudent(p, ws, two_tail)

    def find_plateau(self, x):
        """
        Find a plateau in the middle of a sigmal.
        The plateau is found in the middle of the signal
        :param x:
        :return:
        """
        y = self.pds.find_curve_changes(x)
        y0 = running_mean(y, self.ws)
        runs = zero_runs(y0)
        plat = runs[1]
        i_first = plat[0]
        i_second = plat[1]
        return y0, i_first, i_second

#########################################################################################

def main():
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.style.use('ggplot')
    mpl.rcParams["mathtext.fontset"] = "stix"
    from bessel_filt import bessel_filt_low

    # generate curve
    n = 1024
    x = np.ones(n)
    nx = 200
    start = int(n/2)-int(nx/2)
    end = start + nx
    x[start:end] = 10.0
    fs = 120  # Hz
    dt = 1/fs
    T = dt*n
    t = np.linspace(0, T, n)
    fcut = 8
    order = 5
    p = 0.01    # detect at 5% level (95% confidence)
    ws = 31

    # bessel filter to curve the ends
    xin = bessel_filt_low(x, fcut, fs, order)

    # find the plateau
    fp = FindPlateau(p, ws)
    y0, i_first, i_second = fp.find_plateau(x)

    t1 = i_first*dt
    t2 = i_second*dt

    block = True
    plt.figure()
    plt.plot(t, xin)
    plt.plot(t, y0)
    plt.axvline(x=t1)
    plt.axvline(x=t2)
    plt.show(block)

    # block = True
    # plt.figure()
    # plt.plot(y0)
    # plt.show(block)

    # block = False
    # plt.figure()
    # plt.plot(t, y0)
    # plt.show(block)

    # block = True
    # plt.figure()
    # plt.plot(out)
    # plt.show(block)


if __name__ == '__main__':
    main()









