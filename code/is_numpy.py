import numpy as np
from strcmpi import strcmpi


def is_numpy(a):
    """
    Checks to see if a type is a numpy array or not
    RETURNS
    True = type is a numpy array
    False = type is not a numpy array
    """
    if strcmpi(type(a).__module__, 'numpy') and a.size > 1:
        return True
    else:
        return False


def test_func():
    a = 5
    b = np.asarray([1, 2, 3, 4])
    print(is_numpy(a))
    print(is_numpy(b))


def main():
    test_func()


if __name__ == '__main__':
    main()


