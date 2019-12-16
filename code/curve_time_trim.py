from get_size import length


def curve_time_trim(t, dT):
    """
    For the heat pulse probe, trim the time and dT vectors
    :param t:           as the time vector
    :param dT:          as the change in temperature
    :return:
    """
    t_trim0 = t - t[0]          # start the timestep at t = 0
    t_trim1 = t_trim0[1:]       # trim the first element (cannot use models at time t = 0)
    dT1 = dT[1:]                # trim the first element to be coincident with the start of the timestep
    return t_trim1, dT1         # return the time vector and the change in temperature vector

