

def float_round(x, dec_places):
    """
    Helper function to round a float to a number of decimal places
    """
    thresh = 1.0/(10**dec_places)
    ee = 'f}'
    if x < thresh:
        ee = 'E}'
    fs = '{0:.' + str(int(dec_places)) + ee
    out = float(fs.format(x))
    return out


def float_round_str(x, dec_places):
    """
    Wrapper function for float rounding (return number as string)
    :param x:               as the number
    :param dec_places:      as the number of decimal places
    :return:
    """
    return str(float_round(x, dec_places))
