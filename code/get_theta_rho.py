from constants import *
cw_inv = (1.0/Cw)

"""
HPP experiments to obtain the Liquid Water Content (LWC) and density
"""

def get_theta_rho(k, alpha, theta_o, theta_m, Cm_set=None, Co_set=None):
    """
    Get theta_w and rho
    :param C:               volumetric heat capacity
    :param theta_o:         organic fraction content
    :return:
    """
    C = k / alpha
    theta_w = get_theta_w_from_C_and_organic_content(C, theta_o, theta_m)
    rho = theta_o*rho_o + theta_m*rho_m + theta_w*rho_w
    return theta_w, rho


def constrain_theta(x):
    if x > 1.0:
        return 1.0
    elif x < 0.0:
        return 0.0
    return x


def get_theta_w_from_C_and_organic_content(C, theta_o, theta_m,  Cm_set=None, Co_set=None):
    """
    Obtain the water content from the volumetric heat content and organic content
    :param C:
    :param theta_o:
    :return:
    """
    if Cm_set is not None:
        Cmi = Cm_set
    else:
        Cmi = Cm
    if Co_set is not None:
        Coi = Co_set
    else:
        Coi = Co
    theta_w = cw_inv * (C - (theta_m * Cmi + theta_o * Coi))
    return constrain_theta(theta_w)

