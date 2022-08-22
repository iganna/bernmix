"""
This module includes tools to compute PMFs and CDFs for weighted sum of Bernoulli RVs
"""

import numpy as np

from . import bernmix_control as control
from . import bernmix_int as bmi
from . import bernmix_double as bmd
from . import bernmix_fancy as bm_new


def get_summary(pmf_bm, sum_bias, outcomes = None):
    """
    Prepare output data for PMF
    :param pmf_bm:
    :param sum_bias:
    :param outcomes:
    :return:
    """
    if outcomes is None:
        values = [v + sum_bias for v in range(0, len(pmf_bm))]
        return pmf_bm, values
    else:
        return pmf_bm[outcomes - sum_bias], outcomes


def normalise_params(probs, weights):
    """
    This function normalises parameters of probabilities and weights
    as to
    (1) make weights positive
    (2) remove trivial terms, probability = 0
    (3) remove shifting terms, probability = 1
    :param probs: vector of probabilities
    :param weights: vector of weights
    :return: tuple of * new probabilities
                      * new weights
                      * bias
    """

    # ----------------------------------------------
    # Control Input values
    # ----------------------------------------------
    control.weights_dbl(weights)
    control.probs(probs)
    control.lengths(weights, probs)
    # ----------------------------------------------

    sum_bias = sum([w for w in weights if w < 0])

    # BRV with negative weight is transformed to BRV with opposite probability
    probs = [p if w > 0 else 1 - p
                      for p, w in zip(probs, weights)]
    weights = [abs(w) for w in weights]

    # Remain only significant terms:
    # if some weights or probabilities equal to zero - remove
    idx_significant = [False if (p == 0) | (w == 0) else True for p, w in zip(probs, weights)]


    probs = [p for p, i in zip(probs, idx_significant) if i]
    weights = [w for w, i in zip(weights, idx_significant) if i]

    # Remove constant RVs : BRVs with probabilities equal to 1
    # and change bias
    idx_const = [p == 1 for p in probs]
    sum_bias += sum([w for w, i in zip(weights, idx_const) if i])
    probs = [p for p, i in zip(probs, idx_const) if not i]
    weights = [w for w, i in zip(weights, idx_const) if not i]

    return list(probs), list(weights), sum_bias


# def pmf_int_vals(probs, weights):
#     """
#     This function returns the PMF of the weighted sum of BRVs
#     when weights are integer
#     :param probs:
#     :param weights:
#     :param outcomes:
#     :return: The PMF across all possible values
#     """
#
#     # ----------------------------------------------
#     # Control Input values
#     # ----------------------------------------------
#     control.weights_dbl(weights)
#     control.probs(probs)
#     control.lengths(weights, probs)
#     # ----------------------------------------------
#
#     # remove trivial terms
#     probs, weights, sum_bias = normalise_params(probs, weights)
#     pmf_bm = bmi.pmf(probs, weights)
#
#     values = list(range(0, len(pmf_bm)))
#     values = [v + sum_bias for v in values]
#
#     return pmf_bm, values



def pmf_int(probs, weights, outcomes=None):
    """
    This function returns the PMF of the weighted sum of BRVs
    when weights are integer
    :param probs:
    :param weights:
    :param outcomes:
    :return: The PMF across all possible values
    """

    # ----------------------------------------------
    # Control Input values
    # ----------------------------------------------
    control.weights_dbl(weights)
    control.probs(probs)
    control.lengths(weights, probs)
    # ----------------------------------------------

    # remove trivial terms
    probs, weights, sum_bias = normalise_params(probs, weights)
    pmf_bm = bmi.pmf(probs, weights)

    return get_summary(pmf_bm, sum_bias, outcomes)


def pmf_int_conv(probs, weights, outcomes=None):
    """
    This function returns the PMF of the weighted sum of BRVs
    when weights are integer
    :param probs:
    :param weights:
    :param outcomes:
    :return: The PMF across all possible values
    """

    # ----------------------------------------------
    # Control Input values
    # ----------------------------------------------
    control.weights_dbl(weights)
    control.probs(probs)
    control.lengths(weights, probs)
    # ----------------------------------------------

    # remove trivial terms
    probs, weights, sum_bias = normalise_params(probs, weights)
    pmf_bm = bm_new.dp_int(probs, weights)

    return get_summary(pmf_bm, sum_bias, outcomes)


def pmf_int_conv_fast(probs, weights, outcomes=None):
    """
    This function returns the PMF of the weighted sum of BRVs
    when weights are integer
    :param probs:
    :param weights:
    :param outcomes:
    :return: The PMF across all possible values
    """

    # ----------------------------------------------
    # Control Input values
    # ----------------------------------------------
    control.weights_dbl(weights)
    control.probs(probs)
    control.lengths(weights, probs)
    # ----------------------------------------------

    # remove trivial terms
    probs, weights, sum_bias = normalise_params(probs, weights)
    pmf_bm = bm_new.dp_int_fast(probs, weights)

    return get_summary(pmf_bm, sum_bias, outcomes)


def pmf_int_dep(probs, weights, cov):
    """
    This function returns the PMF of the weighted sum of BRVs
    when weights are integer
    :param probs: pro
    :param weights:
    :param cov: covariance matrix
    :return: The PMF across all possible values
    """

    # ----------------------------------------------
    # Control Input values
    # ----------------------------------------------
    control.weights_dbl(weights)
    control.probs(probs)
    control.lengths(weights, probs)
    # ----------------------------------------------

    # remove trivial terms
    probs, weights, sum_bias = normalise_params(probs, weights)
    pmf_bm = bm_new.dyn_tree_cov(probs, weights, cov)

    return get_summary(pmf_bm, sum_bias)

def cdf_int(probs, weights, outcomes=None):
    """
    This function returns the CDF of the weighted sum of BRVs
    when weights are integer
    when weights are integer
    :param probs:
    :param weights:
    :param outcomes:
    :return: The PMF across all possible values
    """

    # ----------------------------------------------
    # Control Input values
    # ----------------------------------------------
    control.probs(probs)
    control.weights_dbl(weights)
    control.lengths(weights, probs)
    # ----------------------------------------------

    # remove trivial terms
    probs, weights, sum_bias = normalise_params(probs, weights)
    # compute PMF
    pmf_bm = bmi.pmf(probs, weights)
    cdf_bm = np.cumsum(pmf_bm)

    return get_summary(cdf_bm, sum_bias, outcomes)


def cdf_double(probs, weights, target_indiv,
               m_rounding=10**6, n_solutions=None):
    """
    This function reputrn the vector of probabilities for possible values
    of the weighted sum of Bernoulli random variables
    when weights are double
    :param probs: a List of real numbers in the range [0,1]
                 representing probabilities of BRVs
    :param weights: a List of numbers; weights in a weighted sum of BRVs
    :param target_indiv: a List of binary outcomes of BRVs, 0/1 numbers
    :param m_rounding: a number of points to approximate
                       the weighted sum of BRVs
    :param n_solutions: a number of runs for linear integer programming
                        to correct the CDF value
    :return: an approximated or corrected CDF value for the target_indiv
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

    if n_solutions is None:
        cdf_value = bmd.cdf_rounded(probs, weights, target_indiv, m_rounding)
    else:
        cdf_value = bmd.cdf_corrected(probs, weights, target_indiv,
                                  m_rounding, n_solutions)

    return cdf_value


def cdf_permut(probs, weights, target_indiv, n_permut=10 ** 6):
    """
    Get CDF py permutations/simulations
    :param probs: A list of real numbers in the range [0,1]
                 representing probabilities of BRVs
    :param weights: A list of numbers; weights in a weighted sum of BRVs
    :param target_indivs: A list of individual, where each is a list of binary
                          outcomes of BRVs, list with 0 or 1 numbers
    :param n_permut: Number of permutations
    :return: CDF value
    """

    # ----------------------------------------------------
    # Control Input values
    # ----------------------------------------------------
    control.weights_dbl(weights)
    control.probs(probs)
    control.individual(target_indiv)
    control.lengths(weights, probs, target_indiv)
    control.n_solutions(n_permut)
    # ----------------------------------------------------

    # generate a population of Multivariate BRV with target probabilities
    pop = np.zeros((len(probs), n_permut), int)  # Pre-allocate matrix
    for i, p in enumerate(probs):
        pop[i] = np.random.binomial(1, p, n_permut)
    pop = np.transpose(pop)

    # Compute outcomes of the weighted sum of BRVs for the population
    pop_values = list(map(lambda indiv: np.dot(weights, indiv), pop))

    # Compute outcomes of the weighted sum of BRVs for the target individuals
    target_value = np.dot(weights, target_indiv)

    # Compute the approximation of CDF for target_values
    cdf_value = sum(pop_values <= target_value) / n_permut

    return cdf_value


def pmf_int_bf(probs, weights):
    """
    This function caclulates pmf fpr each individual by brute-force search
    :param probs: A list of real numbers in the range [0,1]
                 representing probabilities of BRVs
    :param weights: A list of numbers; weights in a weighted sum of BRVs
    :return: pmf_bm - a vector with discrete values of PMF: P(x) = pmf_bm[x]
    """

    # ----------------------------------------------
    # Control Input values
    # ----------------------------------------------
    control.weights_dbl(weights)
    control.probs(probs)
    control.lengths(weights, probs)
    # ----------------------------------------------

    prob_indiv, outcomes = conv_all_outcomes(probs, weights)
    outcomes_unique = list(set(outcomes))
    outcomes_unique.sort()
    pmf_bm = [sum(prob_indiv[outcomes == i])
              for i in outcomes_unique]
    return pmf_bm, outcomes_unique


def conv_all_outcomes(probs, weights):
    """
    This function calculates probabilities for each outcome
    :param probs: vector of probabilities
    :param weights: vector of weights
    :return: probabilities and outcomes
    """

    def comp_indiv_prob(indiv, probs):
        """ Compute probability for an individual """
        prob_multiply = list(map(lambda ind, p: p if ind == 1 else (1 - p),
                                 indiv, probs))
        return np.prod(prob_multiply)

    # ----------------------------------------------
    # Control Input values
    # ----------------------------------------------
    control.weights_dbl(weights)
    control.probs(probs)
    control.lengths(weights, probs)
    # ----------------------------------------------

    n_terms = len(probs)
    # initialise size of outcomes
    outcomes = np.zeros(2 ** n_terms)
    prob_indiv = np.zeros(2 ** n_terms)

    # initialise two first of outcomes (0,0,...0) and (1,0,...0)
    outcomes[0:2] = [0, weights[0]]
    prob_indiv[0] = comp_indiv_prob(np.zeros(n_terms), probs)
    prob_indiv[1] = comp_indiv_prob(np.append([1], np.zeros(n_terms-1)), probs)

    for i in range(1, n_terms):
        n = 2 ** i
        outcomes[n:2 * n] = outcomes[0:n] + weights[i]
        prob_indiv[n:2 * n] = prob_indiv[0:n] / (1 - probs[i]) * probs[i]

    return prob_indiv, outcomes



def pmf_int_joint(probs, w_vec, v = None, tol = 1e-5):
    """
    Function to predict joint distribution of
    :param p: pribabilities of Bernoulli RVs
    :param w_vec: vector of vectors of weights for RVs
    :param tol: tolerance
    :return:
    """
    # Checks
    # TODO

    # -----------------------------------
    # Preparation
    # -----------------------------------

    # Algorithm


    w_max = [sum(w) + 1 for w in w_vec]
    m_max = np.cumprod(w_max)

    w_new = w_vec[0]
    for i in range(1, len(w_vec)):
        w_new = [w1 + w2 * m_max[i-1] for w1, w2 in zip(w_new, w_vec[i])]

    if v is None:
        s = pmf_int(probs, w_new)
    else:
        s = pmf_int_dep(probs, w_new, v)


    res = [[p, val] for p, val in zip(s[0], s[1]) if p > tol]
    prob_sum = sum([p for p, _ in res])
    res = [[p / prob_sum, val] for p, val in res]

    res_brv = []
    for p, val in res:
        y = []
        for j in range(len(w_max)):
            y_tmp = val % w_max[j]
            y += [y_tmp]
            val -= y_tmp
            val /= w_max[j]
            val = round(val)
        res_brv += [[p] + [y]]

    pmf = [p for p, v in res_brv]
    outcomes = [v for p, v in res_brv]

    return pmf, outcomes


def pmf_int_prod(probs, w_vec, v=None, tol=1e-5):

    pmf_joint, outcomes_joint = pmf_int_joint(probs, w_vec, v, tol)
    pmf_joint = np.array(pmf_joint)
    outcomes_joint = np.array(outcomes_joint)

    outcomes = [np.prod(vals) for vals in outcomes_joint]
    outcomes_unique = list(set(outcomes))
    outcomes_unique.sort()

    pmf_unique = [sum(pmf_joint[outcomes == i]) for i in outcomes_unique]

    return pmf_unique, outcomes_unique


#
# def poibinmix_pmf_int(probs, wights):
#     pass
#
#
# def poibinmix_cdf_double(probs, wights, target_value, n_points = None):
#     pass
#
#
# def binmix_pmf_int(probs, num_of_trails, wights):
#     pass
#
#
# def binmix_cdf_double(probs, num_of_trails, wights, target_value,
#                       n_points = None):
#     pass
#
#
# def radmix_pmf_int(probs, wights):
#     pass
#


if __name__ == "__main__":
    pass
