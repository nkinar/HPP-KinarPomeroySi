from get_size import length


def cut_to_zero(x, t):
    """
    Helper function to cut a sequence to zero
    :param x:       as the sequence
    :param t:       as the associated timestep
    :return:
    """
    n = length(x)
    if n != length(t):
        raise ValueError('cut_to_zero: x and t must have the same length')
    elem = x >= 0
    xout = x[elem]
    tout = t[elem]
    return xout, tout

