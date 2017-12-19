"""

This module generates data for verifying run-time complexity of PMF calculation
and measures the run-time

"""

__author__ = "Anna Igolkina"
__all__ = ["gen_special_int_distr", ]


import glob
import os
from multiprocess import Pool
import numpy as np
import math
import scipy.stats as st
import itertools as it
import bernmix_double.bernmix_int as bmi
from tests.poibin import PoiBin


def error_mn(p1, p2):
    """
    Maximum norm (MN) of difference between distributions
    """
    return max(abs(p1-p2))

def error_rss(p1, p2):
    """
    Root sum of squares (RSS)
    """
    return math.sqrt(sum((p1 - p2) ** 2))


def accuracy_binom_distr(n_range,
                         n_repeats):
    """

    :param n_range:
    :param n_repeats:
    :return:
    """
    for n in n_range:
        mn = 0
        rss = 0
        for _ in range(n_repeats):
            weights = np.ones(n)
            p = np.random.rand(1)
            probs = np.repeat(p, n)
            pmf_exact = binom.pmf(range(0, n + 1), n, p)
            pmf_bermnix = bmi.pmf(probs, weights)
            mn += error_mn(pmf_bermnix, pmf_exact)
            rss += error_rss(pmf_bermnix, pmf_exact)

        yield n, mn / n_repeats, rss / n_repeats


def accuracy_poibin_distr(n_range,
                          n_repeats):
    """

    :param n_range:
    :param n_repeats:
    :return:
    """
    for n in n_range:
        mn = 0
        rss = 0
        for _ in range(n_repeats):
            weights = np.ones(n)
            probs = np.random.rand(n)
            d = PoiBin(probs)
            pmf_exact = d.pmf(range(0, n + 1))
            pmf_bermnix = bmi.pmf(probs, weights)
            mn += error_mn(pmf_bermnix, pmf_exact)
            rss += error_rss(pmf_bermnix, pmf_exact)
        yield n, mn / n_repeats, rss / n_repeats


def gen_special_int_distr(n_range,
                          n_repeats,
                          path_to_out_folder):
    """

    :param n_range:
    :param n_repeats:
    :param path_to_out_folder:
    :return:
    """

    def unif(n):
        """ Uniform distribution """
        return np.random.randint(1, 10, n)

    def even_inif(n):
        """ Even distribution """
        return np.random.randint(1, 5, n) * 2

    def poisson(n):
        """ Poisson distribution """
        w = np.random.poisson(5, n)
        while sum(w == 0) > 0:
            w = np.random.poisson(5, n)
        return w

    def geometr(n):
        """ Geometric Distribution """
        return np.random.geometric(0.2, n) + 1

    def betabin(n):
        """ BetaBinomial Distribution """
        return [st.binom.rvs(9, st.beta.rvs(0.2, 0.2)) + 1 for _ in range(n)]

    # distributions to generate weights
    gen_distr_w = [unif, even_inif, poisson, geometr, betabin]

    n_digits = len(str(n_repeats))
    for n, cnt in it.product(n_range, range(n_repeats)):
        for distr_id, distr_w in enumerate(gen_distr_w):
            probs = np.random.rand(n)
            weights = distr_w(n)
            pw = np.concatenate(([probs], [weights]), axis=0)
            np.savetxt(path_to_out_folder + 'distr_n' +
                       str(n) + '_' +
                       'w' + str(distr_id + 1) + '_' +
                       str(cnt).zfill(n_digits) + '_pw.txt', pw, fmt='%.16f')


def exact_special_int_distr(path_disrt_folder, path_exact_folder):

    def comp_indiv_prob(indiv, probs):
        """ Compute probability for individual """
        prob_multiply = list(map(lambda i: probs[i] if indiv[i] == 1 else (1 - probs[i]),
                                 range(30)))
        return np.prod(prob_multiply)

    def exact_pmf(probs, weights):

        n = len(probs)
        # initialise size of outcomes
        outcomes = np.zeros(2 ** n)
        pmf = np.zeros(2 ** n)

        # initialise two first of outcomes (0,0,...0) and (1,0,...0)
        outcomes[0:2] = [0, weights[0]]
        pmf[0] = comp_indiv_prob(np.zeros(n), probs)
        pmf[1] = comp_indiv_prob(np.append([1], np.zeros(n-1)), probs)

        for i in range(1, n):
            print(i)
            n = 2 ** i
            outcomes[n:2 * n] = outcomes[0:n] + weights[i]
            pmf[n:2 * n] = pmf[0:n] / (1 - probs[i]) * probs[i]
        return pmf, outcomes

    def exact_cdf(file):

        pw = np.loadtxt(file_dir_pw + 'distr' + str(i_distr).zfill(2) + '_pw.txt')
        probs = pw[0]
        weights = pw[1]

        target_indivs = np.loadtxt(file_dir_pw + 'distr' + str(i_distr).zfill(2) + '_pop.txt').astype(int)

        print('Proc', i_distr, 'pmf is computed')

        target_values = list(map(lambda x: np.dot(x, weights), target_indivs))

        cdfs = list(map(lambda x: sum(indiv_probs[indiv_values <= x]), target_values))

        file_dir_res = 'res_examples_n0030/'

        np.savetxt(file_dir_res + 'distr' + str(i_distr).zfill(2) + '_exact.txt',
                   cdfs, fmt='%.16f')

        print('Proc', i_distr, 'is ended')

    with Pool(10) as workers:
        pmap = workers.map
        pmap(cdfs, range(n_disrt))


def accuracy_special_int_distr():

    pass




if __name__ == '__main__':
    # generate distribution
    path_to_out_folder = 'tests/data_test_precision/tmp/'



    #accuracy binomial
    for n in [10, 50, 100, 500, 1000, 5000]: