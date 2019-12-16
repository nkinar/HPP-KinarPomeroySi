from os import listdir
from os.path import isfile, join
from natsort import natsorted
import os


def list_files_in_dir(mypath):
    """
    List all files in directory (not sorted by time)
    :param mypath:
    :return:
    REFERENCE: https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory
    """
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    return onlyfiles


def list_files_in_dir_sorted(mypath):
    """
    Obtain a sorted version of the files in the directory.
    This is sorted according to natural order.  Uses the natsorted library.
    :param mypath:
    :return:
    """
    list_dir = list_files_in_dir(mypath)
    list_dir_sorted = natsorted(list_dir)
    return list_dir_sorted


def list_files_in_dir_time_sorted(dirpath):
    """
    List all of the files in a directory sorted by time
    :param dirpath:
    :return:
    REFERENCE: https://stackoverflow.com/questions/168409/how-do-you-get-a-directory-listing-sorted-by-creation-date-in-python
    """
    a = [s for s in listdir(dirpath)
         if os.path.isfile(os.path.join(dirpath, s))]
    a.sort(key=lambda s: os.path.getmtime(os.path.join(dirpath, s)))
    return a
