import re

def get_all_digits_in_string(s):
    """
    Obtain all digits in a string.
    NOTE that this is returned as a string, not an integer.
    :param s:
    :return:
    REFERENCE: https://stackoverflow.com/questions/11339210/how-to-get-integer-values-from-a-string-in-python
    """
    return (re.findall('\d+', s))