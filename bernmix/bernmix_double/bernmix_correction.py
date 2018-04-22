"""
This module includes tools to correct CDF values for a weighted sum
of Bernoulli RVs using linear integer programming
"""

__author__ = "Anna Igolkina"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Anna Igolkina"
__email__ = "igolkinaanna11@gmail.com"

__all__ = ["cdf_corrected"]

import numpy as np
from cvxopt import glpk, matrix

from . import bernmix_double as bmd
from .. import bernmix_int as bmi
from .. import bernmix_control as control


def cdf_corrected(probs, weights, target_indiv,
                  m_rounding=10**6, n_solutions=100):
    """
    This function computes the corrected value of CDF for a
    :param probs: a List of real numbers in the range [0,1]
                 representing probabilities of BRVs
    :param weights: a List of numbers; weights in a weighted sum of BRVs
    :param target_indiv: a List of binary outcomes of BRVs, 0/1 numbers
    :param m_rounding: a number of points to approximate
                       the weighted sum of BRVs
    :param n_solutions: a number of runs for linear integer programming
                        to correct the CDF value
    :return: an approximated and corrected CDF value for the target_indiv
    """

    # ----------------------------------------------
    # Control Input values
    # ----------------------------------------------
    control.weights_dbl(weights)
    control.probs(probs)
    control.individual(target_indiv)
    control.lengths(weights, probs, target_indiv)
    control.m_rounding(m_rounding, weights)
    control.n_solutions(n_solutions)
    # ----------------------------------------------

    n = len(weights)
    w = weights
    wc = bmd.weight_rounded(w, m_rounding)

    # S and Z values according to the article
    target_value_s = np.dot(w, target_indiv)
    target_value_z = np.dot(wc, target_indiv).astype(int)

    # number of points into two directions
    k = 10
    target_range = np.arange(target_value_z - k, target_value_z + k + 1, 1)
    distr_in_probs = np.zeros((2 * k + 1, 2))
    for i, target_value in enumerate(target_range):
        pop = binprog_multisol(wc, target_value, n_solutions, 10)
        if pop is None:
            distr_in_probs[i, 0] = 0
            distr_in_probs[i, 1] = 0
            continue

        s_values = list(map(lambda x: np.dot(x[0:n], w), pop))

        indiv_probs = list(map(lambda x: comp_indiv_prob(x[0:n], probs), pop))
        sum_prob = sum(indiv_probs)

        distr_in_probs[i, 0] = sum(np.extract(s_values <= target_value_s,
                                              indiv_probs)) / sum_prob
        distr_in_probs[i, 1] = sum(np.extract(s_values > target_value_s,
                                              indiv_probs)) / sum_prob

    pmf = bmi.pmf(probs, wc)

    pmf_range = list(map(lambda t: pmf[t], target_range))
    cdf2 = sum(pmf[:target_value_z + 1 - k - 1:]) + \
           np.dot(pmf_range, distr_in_probs[:, 0])

    return cdf2


def comp_indiv_prob(probs, indiv):
    """ Compute probability for an individual
    :param probs: a List of real numbers in the range [0,1]
                 representing probabilities of BRVs
    :param indiv: a List of binary outcomes of BRVs, 0/1 numbers
    """
    prob_multiply = list(map(lambda i, p: p if i == 1 else (1 - p),
                             indiv, probs))
    return np.prod(prob_multiply)


def binprog_multisol(weights, target_value, n_solutions, n_fixed=10):
    """
    This function returns multiple solutions of
    zero-one linear programming problem

    :param weights: a List of numbers; weights in a weighted sum of BRVs
    :param target_value: an outcome of the weighted sum of BRVs
    :param n_solutions: a number of runs for linear integer programming
                        to correct the CDF value
    :param n_fixed: an area around target_value utilised for correction
    :return: a
    """

    glpk.options['msg_lev'] = 'GLP_MSG_OFF'

    n = len(weights)
    c_init = np.append(weights, -target_value).astype(int)
    c = c_init.reshape((1, n+1))

    a_ident = np.identity(n + 1).astype(int)
    a_ident = np.delete(a_ident, n, axis=0)
    a_ub = np.concatenate((a_ident, -c), axis=0)
    b_ub = np.append(np.ones(n), [0]).astype(int).reshape(n+1, 1)

    pop = np.empty((0, n+1))
    for _ in range(n_solutions):
        # x_(n+1) = 1
        a_eq = np.append(np.zeros(n), [1]).astype(int).reshape(1, n+1)
        b_eq = np.asmatrix([1])

        idx = np.random.permutation(n)
        idx = idx[0:round(n_fixed/100 * n)]
        for i in idx:
            tmp = np.zeros(n+1)
            tmp[i] = 1
            a_eq = np.concatenate((a_eq, [tmp]), axis=0)
            b_eq = np.concatenate((b_eq, [np.random.randint(0, 2, 1)]), axis=0)

        c_mx = matrix(c, tc='d')
        g_mx = matrix(a_ub, tc='d')
        h_mx = matrix(b_ub, tc='d')
        a_mx = matrix(a_eq, tc='d')
        b_mx = matrix(b_eq, tc='d')

        status, x = glpk.ilp(c=c_mx,
                             G=g_mx,
                             h=h_mx,
                             A=a_mx,
                             b=b_mx,
                             B=set(range(n+1)))
        if status != 'optimal':
            continue

        pop = np.concatenate((pop, np.matrix(x).reshape((1, n+1))), axis=0)

    if pop.shape[1] == 0:
        return None

    pop = np.unique(pop, axis=0)
    # ANNA: from 0 to n
    return pop


if __name__ == "__main__":
    print('BernMix_correction is loaded')
