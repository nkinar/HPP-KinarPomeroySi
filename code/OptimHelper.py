"""
Class to minimize a function where the known is used.
This is good if you want to obtain an inverse of a function.
"""


class OptimHelper:

    def __init__(self, f, known):
        """
        Initialize the function
        :param f:
        :param known:
        """
        self.f = f
        self.known = known

    def set_known(self, known):
        """
        Set the known inverse of the function
        :param known:
        :return:
        """
        self.known = known

    def fcall(self, x):
        """
        Function to call when the value is to be calculated.
        :param x:
        :return:
        """
        val = self.f(x) - self.known
        if val < 0:     # take absolute value
            out = -val
        else:
            out = val
        return out

