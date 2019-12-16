#!/usr/bin/env python3


def is_even(x):
    return isEven(x)


def is_odd(x):
    return isOdd(x)


def isEven(x):
    if x % 2 == 0:
        return True
    else:
        return False


def isOdd(x):
    return not isEven(x)


# TEST FUNCTIONS
def test_even_odd():
    x = 5
    y = 6
    print('x = ', x, isOdd(x))
    print('y = ', y, isEven(y))
    print('0 is even:', isEven(0))
    

def main():
    test_even_odd()
    
if __name__ == '__main__': main()