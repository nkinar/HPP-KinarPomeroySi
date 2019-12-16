#!/usr/bin/env python3
import numpy as np

# REFERENCES:
# [1] Knuth D. The Art of Computer Programming, vol II.
# [2] https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison/17341


def get_machine_epsilon():
    return np.finfo(float).eps

def get_float32_epsilon():
    return np.finfo(np.float32).eps

EPS = get_machine_epsilon()

def approximately_equal(a, b, epsilon=EPS):
    lhs = np.abs(a-b)
    abs_a = np.abs(a)
    abs_b = np.abs(b)
    if abs_a < abs_b:  # less than <
        rhs = abs_b
    else:
        rhs = abs_a
    rhs *= epsilon
    return lhs <= rhs


def essentially_equal(a, b, epsilon=EPS):
    lhs = np.abs(a-b)
    abs_a = np.abs(a)
    abs_b = np.abs(b)
    if abs_a > abs_b:  # greater than >
        rhs = abs_b
    else:
        rhs = abs_a
    rhs *= epsilon
    return lhs <= rhs


def definitely_greater_than(a, b, epsilon=EPS):
    lhs = a-b
    abs_a = np.abs(a)
    abs_b = np.abs(b)
    if abs_a < abs_b:
        rhs = abs_b
    else:
        rhs = abs_a
    rhs *= epsilon
    return lhs <= rhs


def definitely_less_than(a, b, epsilon=EPS):
    lhs = b-a
    abs_a = np.abs(a)
    abs_b = np.abs(b)
    if abs_a < abs_b:
        rhs = abs_b
    else:
        rhs = abs_a
    rhs *= EPS
    return lhs <= rhs


def almost_equal(a, b, epsilon=EPS):
    return essentially_equal(a, b, epsilon)


def very_close(a, b, epsilon=EPS):
    approximately_equal(a, b, epsilson)


def close(a, b, epsilon=EPS):
    return approximately_equal(a, b, epsilon)