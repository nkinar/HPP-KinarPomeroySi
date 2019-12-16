from collections import OrderedDict

"""
How to create a 3D dictionary.  For default, this is an ordered dictionary.
REFERENCE: https://stackoverflow.com/questions/12167192/pythonic-way-to-create-3d-dict
"""

class AutoDict(OrderedDict):
    def __missing__(self, key):
        x = AutoDict()
        self[key] = x
        return x