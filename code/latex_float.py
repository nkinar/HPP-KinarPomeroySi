from float_comparison import close

def latex_float(float_str):
    """
    Reformats number as latex float
    REFERENCE: https://stackoverflow.com/questions/13490292/format-number-using-latex-notation-in-python/13490601
    :param f:
    :return:
    """
    if "e" in float_str:
        base, exponent = float_str.split("e")
        exponent = int(exponent)
        if exponent == 0:
            out = r"{0}".format(base)  # reformat as single number
        else:
            # reformat as scientific number
            out = r"{0} \times 10^{{{1}}}".format(base, int(exponent))
        return out
    else:
        return float_str