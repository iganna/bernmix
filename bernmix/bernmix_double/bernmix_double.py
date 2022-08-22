"""
This module contain tools to compute CDF for a weighted sum of BRVs
with non-integer weights using rounding
"""

__author__ = "Anna Igolkina"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Anna Igolkina"
__email__ = "igolkinaanna11@gmail.com"

__all__ = ["weight_rounded", "cdf_rounded"]

import numpy as np

from .. import bernmix_int as bmi
from .. import bernmix_control as control


def weight_rounded(weights, m_points):
    """
    This function return the scaled and rounded  approximation of weights
    :param weights: a List of numbers; weights in a weighted sum of BRVs
    :param m_points: a number of points to approximate
                       the weighted sum of BRVs
    :return:
    """

    # ----------------------------------------------
    # Control Input values
    # ----------------------------------------------
    control.weights_dbl(weights)
    control.m_rounding(m_points, weights)
    # ----------------------------------------------

    n = len(weights)
    weight_sum = sum(weights)
    c = (m_points + n) / weight_sum
    return np.around(weights * c).astype(int)


def cdf_rounded(probs, weights, indiv, m_rounding=10**6):
    """
    This function calculates the CDF value for several individuals
    using fixed number of points to approximate
    :param probs: a List of real numbers in the range [0,1]
                  representing probabilities of BRVs
    :param weights: a List of numbers; weights in a weighted sum of BRVs
    :param target_indivs: a list of individuals where each is a list of binary
                          outcomes of BRVs, 0/1 numbers
    :param m_rounding: a number of points to approximate
                       the weighted sum of BRVs
    :return:
    """

    # ----------------------------------------------
    # Control Input values
    # ----------------------------------------------
    control.weights_dbl(weights)
    control.probs(probs)
    control.individual(indiv)
    control.lengths(weights, probs, indiv)
    control.m_rounding(m_rounding, weights)
    # ----------------------------------------------

    # Compute approximated PMF
    weights_new = weight_rounded(weights, m_rounding)
    pmf = bmi.pmf(probs, weights_new)

    value = np.dot(weights_new, indiv).astype(int)
    cdf_value = sum(pmf[:(value+1)])

    return cdf_value
