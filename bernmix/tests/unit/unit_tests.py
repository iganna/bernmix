"""
Unit tests
"""

__author__ = "Anna Igolkina"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Anna Igolkina"
__email__ = "igolkinaanna11@gmail.com"

import unittest
import random

import numpy as np
from ...bernmix import bernmix_cdf_int, bernmix_pmf_int, bernmix_cdf_double


class TestBernMix(unittest.TestCase):
    def test_int_case(self):

        # define number of terms
        n = random.randint(1, 100)
        weight_sum = random.randint(n, 100)

        # randomly define probabilities of i.ni.d. BRVs
        probs = np.random.rand(n)

        # define weights
        weights = np.random.rand(n)
        weights = np.round(weights * weight_sum / sum(weights)).astype(int)

        # no weights should be equal to zero
        weights[weights == 0] += 1

        # run
        bernmix_pmf_int(probs, weights)
        bernmix_cdf_int(probs, weights)

    def test_double_case(self):

        # define number of terms
        n = random.randint(1, 100)
        weight_sum = random.randint(n, 100)

        # randomly define probabilities of i.ni.d. BRVs
        probs = np.random.rand(n)

        # define weights
        weights = np.random.rand(n)
        weights = np.round(weights * weight_sum / sum(weights)).astype(int)

        # no weights should be equal to zero
        weights[weights == 0] += 1

        # run
        bernmix_cdf_double(probs, weights)


if __name__ == "__main__":
    unittest.main()