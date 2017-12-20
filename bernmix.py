


__author__ = "Anna Igolkina"
__license__ = "MIT"
__version__ = ""
__maintainer__ = "Anna Igolkina"
__email__ = "igolkinaanna11@gmail.com"
__status__ = "Development"

__all__ = ["bernmix_pmf_int", "bernmix_cdf_double"]

import bernmix_double.bernmix_double as bmd
import bernmix_int.bernmix_int as bmi
import numpy as np

def normalise_params(probs, weights):
    """
    This function normalises parameters of probabilities and weights as to make weights positive
    :param probs: vector of probabilities
    :param weights: vector of weights
    :return: sum of negative weights
    """
    sum_neg_weights = sum([w for w in weights if w < 0])

    probs = [p * (w > 0) + (1 - p) * (w < 0) for p, w in zip(probs, weights)]
    weights = [abs(w) for w in weights]

    # if some weights == 0
    # if some probabilities == 0

    return probs, weights, sum_neg_weights

def exact_int_pmf(probs, weights):
    """
    This function caclulates pmf fpr each individual by convolution
    :param probs: vector of probabilities
    :param weights: vector of weights
    :return: probabilities and outcomes
    """
    if len(probs) != len(weights):
        raise ValueError('Vectors of probabilities and weights must be the same length')
    if len(probs) == 0:
        raise ValueError('Sum of BRVs should contain at least one RV')

    prob_indiv, outcomes = exact_prob_indiv(probs, weights)
    pmf_bm = [sum(prob_indiv[outcomes == i]) for i in range(0, sum(weights) + 1)]
    return pmf_bm


def exact_prob_indiv(probs, weights):
    """
    This function caclulates probabilities for each outcome
    :param probs: vector of probabilities
    :param weights: vector of weights
    :return: probabilities and outcomes
    """
    def comp_indiv_prob(indiv, probs):
        """ Compute probability for an individual """
        n = len(probs)
        prob_multiply = list(map(lambda i: probs[i] if indiv[i] == 1 else (1 - probs[i]),
                                 range(n)))
        return np.prod(prob_multiply)

    if len(probs) != len(weights):
        raise ValueError('Vectors of probabilities and weights must be the same length')
    if len(probs) == 0:
        raise ValueError('Sum of BRVs should contain at least one RV')

    n = len(probs)
    # initialise size of outcomes
    outcomes = np.zeros(2 ** n)
    prob_indiv = np.zeros(2 ** n)

    # initialise two first of outcomes (0,0,...0) and (1,0,...0)
    outcomes[0:2] = [0, weights[0]]
    prob_indiv[0] = comp_indiv_prob(np.zeros(n), probs)
    prob_indiv[1] = comp_indiv_prob(np.append([1], np.zeros(n-1)), probs)

    for i in range(1, n):
        n = 2 ** i
        outcomes[n:2 * n] = outcomes[0:n] + weights[i]
        prob_indiv[n:2 * n] = prob_indiv[0:n] / (1 - probs[i]) * probs[i]

    return prob_indiv, outcomes


def bernmix_pmf_int(probs, weights, outcomes=None):
    """
    This function reputrn the vector of probabilities for possible values
    of the weighted sum of Bernoulli random variables
    when weights are integer
    :param probs:
    :param weights:
    :param outcomes:
    :return: The PMF across all possible values
    """
    if len(probs) != len(weights):
        raise ValueError('Vectors of probabilities and weights must be the same length')
    if len(probs) == 0:
        raise ValueError('Sum of BRVs should contain at least one RV')
    if sum((0 > probs) | (probs > 1)):
        raise ValueError('Probabilities should be within [0 1] segment')

    # ANNA
    # if weighas are integer
    # if outcomes are integer

    probs, weights, sum_neg_weights = normalise_params(probs, weights)
    pmf_bm = bmi.pmf(probs, weights)

    if outcomes in None:
        return pmf_bm
    else:
        return pmf_bm[outcomes + sum_neg_weights]


def bernmix_cdf_double(probs, weights, target_indiv, m=10**6, l=100):
    """
    This function reputrn the vector of probabilities for possible values
    of the weighted sum of Bernoulli random variables
    when weights are double
    :param probs:
    :param weights:
    :param target_value:
    :param n_points: a number of point to approximate distribution
    :return: the CDF in target_value
    """
    pass


def poibinmix_pmf_int(probs, wights):
    """
    Linear combination of Poisson Binomial distributions
    :param probs:
    :param wights:
    :return:
    """
    pass


def poibinmix_cdf_double(probs, wights, target_value, n_points = None):
    """
    Linear combination of Poisson Binomial distributions
    :param probs:
    :param wights:
    :param target_value:
    :param n_points:
    :return:
    """
    pass


def binmix_pmf_int(probs, num_of_trails, wights):
    """
    Linear combination of Poisson Binomial distributions
    :param probs:
    :param num_of_trails:
    :param wights:
    :return:
    """
    pass

def binmix_cdf_double(probs, num_of_trails, wights, target_value, n_points = None):
    """
    Linear combination of Poisson Binomial distributions
    :param probs:
    :param num_of_trails:
    :param wights:
    :param target_value:
    :param n_points:
    :return:
    """
    pass

def radmix_pmf_int(probs, wights):
    """
    Rademacher distribution
    :return:
    """
    pass




if __name__ == "__main__":
    pass