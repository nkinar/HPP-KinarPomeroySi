from constants import *


def get_volumetric_heat_capacity(theta_m, theta_o, theta_w):
    """
    Get the volumetric heat capacity of a soil
    :param theta_m:         mineral fraction
    :param theta_o:         organic fraction
    :param theta_w:         water fraction
    Constants are from Van Wijk and DeVries (1963)
    :return:
    """
    C = theta_m*Cm + theta_o*Co + theta_w*Cw
    return C



