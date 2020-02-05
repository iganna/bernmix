"""
This module includes tools to compute PMFs and CDFs for weighted sum of Bernoulli RVs
"""

import numpy as np

from . import bernmix_control as control
from . import bernmix_int as bmi
from . import bernmix_double as bmd


def normalise_params(probs, weights):
    """
    This function normalises parameters of probabilities and weights
    as to make weights positive
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
    probs = np.array([p * (w > 0) + (1 - p) * (w < 0)
                      for p, w in zip(probs, weights)])
    weights = np.array([abs(w) for w in weights])

    # Remain only significant terms:
    # if some weights or probabilities equal to zero - remove
    idx_significant = ~ ((probs == 0) | (weights == 0))
    probs = probs[idx_significant]
    weights = weights[idx_significant]

    # Remove constant RVs : BRVs with probabilities equal to 1
    # and change bias
    idx_const = probs == 1
    sum_bias += sum(weights[idx_const])
    probs = probs[~ idx_const]
    weights = weights[~ idx_const]

    return probs, weights, sum_bias


def pmf_int_vals(probs, weights):
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

    values = list(range(0, len(pmf_bm)))
    values = [v + sum_bias for v in values]

    return pmf_bm, values



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

    if outcomes is None:
        return pmf_bm
    else:
        return pmf_bm[outcomes - sum_bias]


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

    if outcomes is None:
        return cdf_bm
    else:
        return cdf_bm[outcomes + sum_bias]


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


def cdf_permut(probs, weights, target_indivs, n_permut=10 ** 6):
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
    [control.individual(indiv) for indiv in target_indivs]
    control.lengths(weights, probs, *target_indivs)
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
    target_values = map(lambda indiv: np.dot(weights, indiv), target_indivs)

    # Compute the approximation of CDF for target_values
    cdfs = [map(lambda value: sum(pop_values <= value) / n_permut,
                target_values)]

    return cdfs


def conv_pmf_int(probs, weights):
    """
    This function caclulates pmf fpr each individual by convolution
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
    pmf_bm = [sum(prob_indiv[outcomes == i])
              for i in range(0, sum(weights) + 1)]
    return pmf_bm


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
