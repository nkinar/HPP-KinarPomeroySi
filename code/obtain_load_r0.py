import numpy as np
from get_size import length
from check_same_size import check_same_size
from constants import *


class SaveLoadR0:
    ERR_STR = 'SaveLoadR0: q_range and r0 not the same size'

    def __init__(self):
        self.lookup = {}

    def get_r0(self, qval):
        """
        Obtain the r0 value, given the q
        :param qval:
        :return:
        """
        return self.lookup[str(qval)]

    def save_r0(self, path, q_range, r0):
        """
        Save r0 to file
        :param path:
        :param q_range:
        :param r0:
        :return:
        """
        check_same_size(q_range, r0, self.ERR_STR)
        np.savez(path, q_range=q_range, r0=r0)

    def load_r0(self, path):
        """
        Load the radius r0
        :return:
        """
        array = np.load(path)
        q_range = array[QRANGE_IDENT]
        r0 = array[R0_IDENT]
        check_same_size(q_range, r0, self.ERR_STR)
        n = length(q_range)
        for k in range(n):
            s = str(q_range[k])
            self.lookup[s] = r0[k]
    # DONE

