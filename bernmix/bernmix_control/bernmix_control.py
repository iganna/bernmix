"""
This module contain methods to control arguments of functions
"""

__author__ = "Anna Igolkina"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Anna Igolkina"
__email__ = "igolkinaanna11@gmail.com"

__all__ = ["weights_dbl",
           "weights_int",
           "probs",
           "individual",
           "lengths",
           "m_rounding",
           "n_solutions",
           "outcomes_sum"]

import numbers
import collections


def weights_dbl(weights):
    """
    Control whether weights parameter is a List of numbers
    """
    if not isinstance(weights, collections.Iterable):
        raise ValueError('Weights are not Iterable')
    if not all([isinstance(w, numbers.Number) for w in weights]):
        raise ValueError('Weights are not numbers')


def weights_int(weights):
    """
    Control whether weights parameter is a List of integer numbers
    """
    weights_dbl(weights)
    if not all([isinstance(w, int) for w in weights]):
        raise ValueError('Weights are not integer')


def probs(probs):
    """
    Control whether probs parameter is a List of real numbers in [0, 1]
    """
    if not isinstance(probs, collections.Iterable):
        raise ValueError('Probabilities are not Iterable')
    if not all([isinstance(p, numbers.Number) for p in probs]):
        raise ValueError('Probabilities are not numbers')
    if not all([0 <= p <= 1 for p in probs]):
        raise ValueError('Probabilities are not in [0, 1]')


def individual(indiv):
    """
    Control whether an individual is a list of 0 or 1
    """
    if not isinstance(indiv, collections.Iterable):
        raise ValueError('Individual is not Iterable')
    if not all([isinstance(t, int) for t in indiv]):
        raise ValueError('Target individual does not contain all numbers')
    if not all([t in (0, 1) for t in indiv]):
        raise ValueError('Target individual is a vector of [0 and 1]')


def lengths(*args):
    """
    Control whether arguments are lists of the same length
    """
    n = [len(element) for element in args]
    if len(n) == 0:
        raise ValueError('Please provide arguments')
    if not all([n[0] == elem for elem in n]):
        raise ValueError('Parameters should be the same length')


def m_rounding(m, weights):
    """
    Control the M-value for rounding
    """
    if not isinstance(m, int):
        raise ValueError('M for rounding should be an integer number')
    if m < len(weights):
        raise ValueError('M is lower than sum of weights')


def n_solutions(n):
    """
    Control the L-value for correction
    """
    if n is None:
        return
    if not isinstance(n, int):
        raise ValueError('N for linear programming should be an integer number')
    if n <= 0:
        raise ValueError('N for linear programming should be positive')


def outcomes_sum(outcomes, weights):
    """
    Control possible outcomes
    """
    error_message = 'Outcomes of a weighted sum are not correct'
    if isinstance(outcomes, collections.Iterable):
        if not all([isinstance(outcome, numbers.Number) for outcome in outcomes]):
            raise ValueError(error_message)
        if not all([0 <= outcome <= sum(weights) for outcome in outcomes]):
            raise ValueError(error_message)
    elif isinstance(outcomes, numbers.Number):
        if not (0 <= outcomes <= sum(weights)):
            raise ValueError(error_message)
    else:
        raise ValueError(error_message)
