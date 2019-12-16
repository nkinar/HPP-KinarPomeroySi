import collections


def check_if_iterable(obj):
    """
    Function to check if an object is iterable
    :param obj:
    :return:
    """
    try:
        iterator = iter(obj)
    except TypeError:
        return False
    return True


def iterable(arg):
    """
    REFERENCE:
    https://stackoverflow.com/questions/1055360/how-to-tell-a-variable-is-iterable-but-not-a-string
    :param arg:
    :return:
    """
    return (
        isinstance(arg, collections.Iterable)
        and not isinstance(arg, str)
    )
