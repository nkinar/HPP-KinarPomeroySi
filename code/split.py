import re

def split(delimiters, string, maxsplit=0):
    """
    Split a string with multiple delimiters
    :param delimiters:
    :param string:
    :param maxsplit:
    :return:
    REFERENCE: https://stackoverflow.com/questions/4998629/python-split-string-with-multiple-delimiters
    """
    regexPattern = '|'.join(map(re.escape, delimiters))
    return re.split(regexPattern, string, maxsplit)