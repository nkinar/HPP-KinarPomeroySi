

def get_alpha_k_C(k, C):
    """
    Obtain the alpha from the thermal conductivity and volumetric heat capacity
    :param k:   thermal conductivity
    :param C:   volumetric heat capacity
    :return:
    """
    alpha = k / C
    return alpha


def get_alpha(k, rho, c):
    """
    Obtain the alpha from density and heat capacity
    :param k:
    :param rho:
    :param c:
    :return:
    """
    alpha = k / (rho * c)
    return alpha