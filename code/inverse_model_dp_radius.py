import numpy as np
from get_size import length
import sys
from numpy import pi
from scipy.optimize import minimize
from sp_model_signal import sp_model_signal_inv
from butterworth_low import butter_lowpass_filter

sys.path.append('../expi-library')
sys.path.append('../derivative-solver')
from Expi import Expi
exp_integral = Expi()

from deriv_forward_inv import forward_first_derivative, inverse_first_derivative, DerivSolver


def expi_inv_vect_neg(known):
    """
    Take the inverse of expi over a vector.  Here it is known that
    x < 0.

    NOTE that this function cannot be applied if x > 0.

    :param known:
    :param neg:
    :param radius:
    :return:
    """
    n = length(known)
    out = np.zeros(n)
    for k in range(n):
        out[k] = exp_integral.Eid_inv_negative(known[k])
    return out


def optim_func_dual_probe_alpha(x, *coeff):
    """
    Dual probe optimization function to obtain alpha.
    This is only applied after r(t) is known.

    :param x:
    :param coeff:
    :return:
    """
    # coefficients that are passed in (q, k, r, t, gamma3)
    q = coeff[0]                # scalar (heat input)
    k = coeff[1]                # scalar (thermal conductivity)
    r = coeff[2]                # vector (change in radius)
    t = coeff[3]                # time vector (s)
    gamma3 = -coeff[4]          # gamma3 = -(r(t)**2)/(4*alpha*t), but we use the positive version
    # coeff to be found
    alpha = x[0]                # alpha as thermal diffusivity
    # computed temperature change for dual probe
    computed = (r**2)/(4*alpha*t)
    rv = np.sum((gamma3 - computed)**2)
    return rv


def obtain_alpha(q, kfound, r, t, gamma3):
    """
    Obtain thermal diffusivity value.
    This is called only after the radius change is found.

    :param q:
    :param kfound:
    :param gamma9:
    :param t:
    :param gamma3:
    :return:
    """
    coeff = (q, kfound, r, t, gamma3)  # known coefficients
    x0 = np.asarray([1.0])
    res = minimize(optim_func_dual_probe_alpha, x0, args=coeff, method='Nelder-Mead', tol=1.0e-15)
    alpha_found = res.x
    return alpha_found[0]


class InverseWarmupCurveDP:
    def __init__(self):
        self.ds = None

    def _set_ds(self, n):
        """
        Check to see if the derivative solver object needs to be re-created
        :param n:
        :return:
        """
        if self.ds is None:
            self.ds = DerivSolver(n)
        else:
            ncheck = self.ds.get_n()
            if n != ncheck:
                self.ds = DerivSolver(n)

    def inverse_warmup_curve_dp(self, k, dT_dual, q, t, r0, typ, lowpass=False, epsilon=0.0):
        """
        Obtain the inverse of the curve data
        Note that the inverse solution cannot be used for t = 0 since the solution is undefined at this point
        The first element is therefore sampled at timestep t = delta_t

        :param k:           as the thermal conductivity determined from the SP model (or assumed)
        :param dT_dual:     dual-probe warmup curve at distance of r(t) (Celcius or K)
        :param q:           inferred heat input (W/m)
        :param t:           vector of regular timestamps (t > 0) for DP
        :param r0:          initial single-probe radius

        :param typ:         'iterative' to use an iterative solution
                            'matrix' to use a matrix solution
                            'regularized' to use Tikhonov regularization

        :param epsilon:     regularization parameter to use in the Tikhanov solution

        :return: (alpha, r_t)
        alpha = thermal diffusivity (m^2 / s)
        r_t = time series showing change in r with respect to time t
        """
        n = length(t)
        dt = t[1] - t[0]    # timestep must be regular
        gamma2 = -dT_dual * ((4 * pi * k) / q)
        gamma3 = expi_inv_vect_neg(gamma2)
        gamma4 = -gamma3
        gamma5 = np.sqrt(gamma4 * t)
        gamma6 = np.log(gamma5)
        gamma7 = forward_first_derivative(gamma6, dt)
        log_r = np.log(r0)
        if typ == 'iterative':
            gamma8 = inverse_first_derivative(gamma7, dt, log_r)
        elif typ == 'matrix':
            self._set_ds(n)
            gamma8 = self.ds.solve_system_normal(gamma7, dt, log_r)
        elif typ == 'regularized':
            self._set_ds(n)
            gamma8 = self.ds.solve_system_tikreg_sparse(gamma7, dt, log_r, epsilon)
        else:
            raise ValueError('inverse_warmup_curve: the type of inverse has not been specified')
        gamma9 = np.exp(gamma8)
        r_t = gamma9
        alpha = obtain_alpha(q, k, r_t, t, gamma3)
        return alpha, r_t


def main():
    pass


if __name__ == '__main__':
    main()

