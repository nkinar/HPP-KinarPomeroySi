from numpy import genfromtxt
from constants import *


def construct_load_file_string(path, start_string, q, t_heat, t_total, additional_text):
    """
    Construct the file name from which to load the data
    :param path:                as the path to the data directory
    :param start_string:        starting part of the file name (different for calibration file)
    :param q:                   heat input [s]
    :param t_heat:              time of heating
    :param t_total:             total time of the experiment
    :param additional_text:     additional text at the end of the file name
    :return:
    NOTE that the heat input is given in seconds, not in minutes.
    """
    fn = path + start_string + str(int(q)) + WM_STR + str(int(t_heat)) + SEC_STR + HEATING_STR + str(int(t_total)) \
         + MIN_TOTAL_STR + additional_text + CSV_EXT
    return fn


def load_data_file(full_path):
    """
    Load a data file given the full path to the file
    :param full_path:           as the full path to the file
    :return:
    """
    data = genfromtxt(full_path, dtype='str', delimiter=DATA_DELIMITER)
    num = data[:, NUM_COL].astype(int)
    Vout = data[:, VOUT_COL].astype(float)
    dV = data[:, dV_COL].astype(float)
    I = data[:, I_COL].astype(float)
    Rw = data[:, RW_COL].astype(float)
    T1 = data[:, T1_COL].astype(float)
    T2 = data[:, T2_COL].astype(float)
    dac_code = data[:, DAC_CODE_COL].astype(int)
    qprime = data[:, QPRIME_CODE_COL].astype(float)
    return num, Vout, dV, I, Rw, T1, T2, dac_code, qprime



