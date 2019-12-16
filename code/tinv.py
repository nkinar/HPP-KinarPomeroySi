from scipy import stats

"""
Student's t-test significance values

NOTES: 
1. The significance level is given as a frational percentage (0.05 for 5%)
This is similar to the Ebdon book.

REFERENCE:
https://stackoverflow.com/questions/19339305/python-function-to-get-the-t-statistic
"""


def tinv_two_tail(p, df):
    """
    Compute the two tail distribution for Student's t
    :param p:       as the significance level
    :param df:      as the degrees of freedom
    :return:
    """
    pp = 1.0 - 0.5*p
    return stats.t.ppf(pp, df)


def tinv_one_tail(p, df):
    """
    Compute the one tail distribution for Student's t
    :param p:       as the significance level
    :param df:      as the degrees of freedom
    :return:
    """
    pp = 1.0 - p
    return stats.t.ppf(pp, df)


def main():
    p = 0.05
    df = 13
    out = tinv_two_tail(p, df)
    print(out)


if __name__ == '__main__':
    main()