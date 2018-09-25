"""

This module generates data for verifying run-time complexity of PMF calculation
and measures the run-time

"""

__author__ = "Anna Igolkina"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Anna Igolkina"
__email__ = "igolkinaanna11@gmail.com"

__all__ = ["accuracy_binom_distr",
           "accuracy_poibin_distr",
           "gen_special_int_distr",
           "accuracy_special_int_distr"]


import re
import glob
import os

import numpy as np
import math
import scipy.stats as st
import itertools as it
from multiprocess import Pool

from ...bernmix import bernmix_cdf_int, bernmix_pmf_int
from .poibin import PoiBin
from bernmix import conv_pmf_int


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
    This function calculate MN and RSS discrepancy between PMFs for the Binomial distribution
    calculated via theoretical formula and bernmix package
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
            pmf_exact = st.binom.pmf(range(0, n + 1), n, p)
            pmf_bermnix = bernmix_pmf_int(probs, weights)
            mn += error_mn(pmf_bermnix, pmf_exact)
            rss += error_rss(pmf_bermnix, pmf_exact)

        yield '{}\t{:.16f}\t{:.16f}'.format(n, mn / n_repeats, rss / n_repeats)


def accuracy_poibin_distr(n_range,
                          n_repeats):
    """
    This function calculate MN and RSS discrepancy between PMFs
    for the Poisson Binomial distribution
    calculated via theoretical formula and bernmix package
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
            pmf_bermnix = bernmix_pmf_int(probs, weights)
            mn += error_mn(pmf_bermnix, pmf_exact)
            rss += error_rss(pmf_bermnix, pmf_exact)

        yield '{}\t{:.16f}\t{:.16f}'.format(n, mn / n_repeats, rss / n_repeats)


def gen_special_int_distr(n_range,
                          n_repeats,
                          path_distr_folder):
    """
    Generator of distributions, when weights follow one of five distributions:
    Uniform, Even uniform, Poisson, Geometric and Beta-Binomial
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
        for distr_id, distr_func in enumerate(gen_distr_w):
            probs = np.random.rand(n)
            weights = distr_func(n)
            pw = np.concatenate(([probs], [weights]), axis=0)
            np.savetxt(path_distr_folder + 'distr_n' +
                       str(n) + '_' +
                       'w' + str(distr_id + 1) + '_' +
                       str(cnt).zfill(n_digits) + '_pw.txt', pw, fmt='%.16f')



def accuracy_special_int_distr(n_range,
                               path_distr_folder,
                               path_accuracy_folder,
                               n_threads):
    """
    Comparisson of a PMF computed by convolution and computed by bernmix
    :param n_range:
    :param path_distr_folder:
    :param path_accuracy_folder:
    :param n_threads:
    :return:
    """


    def comp_errors(file_pw):
        """ Compute MN and RSS errors for each distribution"""

        # get id of the distribution and id of weights
        n_id, w_id = re.findall(r'/distr_n(\d+)_w(\d+)\w+.txt', file_pw)[0]
        n_id /= 10

        pw = np.loadtxt(file_pw)
        probs = pw[0]
        weights = pw[1]
        pmf_exact = conv_pmf_int(probs, weights)
        pmf_bermnix = bernmix_pmf_int(probs, weights)

        mn = error_mn(pmf_bermnix, pmf_exact)
        rss = error_rss(pmf_bermnix, pmf_exact)
        mn_values[n_id-1, w_id-1] += mn
        rss_values[n_id - 1, w_id - 1] += rss
        comp_repeats[n_id - 1, w_id - 1] += 1


    # initialisation
    n_distributions_for_weights = 5
    mn_values = np.zeros(len(n_range), n_distributions_for_weights)
    rss_values = np.zeros(len(n_range), n_distributions_for_weights)
    comp_repeats = np.zeros(len(n_range), n_distributions_for_weights)


    files = glob.glob(path_distr_folder + '*_pw.txt')
    with Pool(n_threads) as workers:
        pmap = workers.map
        pmap(comp_errors, files)

    mn_values /= comp_repeats
    rss_values /= comp_repeats

    np.savetxt(path_accuracy_folder + 'accuracy_mn.txt', mn_values, '%.16f')
    np.savetxt(path_accuracy_folder + 'accuracy_rss.txt', rss_values, '%.16f')


if __name__ == '__main__':
    pass
    # # Accuracy for Binomial and Poisson Binomial distributions
    # n_range = [10, 50, 100, 500, 1000, 5000]
    # n_repeats = 10
    # results_binom = accuracy_binom_distr(n_range, n_repeats)
    # results_poibin = accuracy_poibin_distr(n_range, n_repeats)
    # # Save results
    # path_accuracy_binpoibin = 'tests/data_test_precision/accuracy_bin_poibin/'
    # if not os.path.exists(path_accuracy_binpoibin):
    #     os.makedirs(path_accuracy_binpoibin)
    # np.savetxt(path_accuracy_binpoibin + 'results_binom.txt', results_binom, '%s')
    # np.savetxt(path_accuracy_binpoibin + 'results_poibin.txt', results_poibin, '%s')
    #
    # # Generate distributions with specifically distributed weights
    # path_distr_folder = 'tests/data_test_precision/distr_int_weights/'
    # n_range = [10, 20, 30]
    # n_repeats = 20
    # gen_special_int_distr(n_range,
    #                       n_repeats,
    #                       path_distr_folder)
    # # Calculate accuracy
    # n_threads = 10
    # accuracy_special_int_distr(n_range,
    #                            path_distr_folder, path_distr_folder,
    #                            n_threads)


