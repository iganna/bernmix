"""

This module generates data for verifying run-time complexity of PMF calculation
and measures the run-time

"""

__author__ = "Anna Igolkina"
__all__ = ["gen_double_distr", ]

import time
import glob
import os
from multiprocess import Pool
import numpy as np
import itertools as it
import bernmix_double.bernmix_int as bmi
from bernmix import exact_pmf


def gen_double_distr(n, n_distr, pop_size,
                     path_to_out_folder):
    """
    This function generates parameters of distributions of positively-weighted sum of BRVs
    with non-integer weights and a population of multivariate-Bernoulli individuals
    :param n_distr:
    :param pop_size:
    :param path_to_out_folder:
    :return:
    """
    n_digits = len(str(n_distr))
    for cnt in range(n_distr):
        pw = np.random.rand(2, n)
        pop = np.random.randint(0, 2, (pop_size, n))
        np.savetxt(path_to_out_folder + 'distr' + str(cnt).zfill(n_digits) + '_pw.txt', pw, fmt='%.16f')
        np.savetxt(path_to_out_folder + 'distr' + str(cnt).zfill(n_digits) + '_pop.txt', pop, fmt='%i')


def exact_cdf_pop(path_disrt_folder, path_cdfs_folder, n_threads=1):
    """
    This function calculate CDF values for several Multivariate Bernoulli RVs
    Paralleled
    :param path_disrt_folder:
    :param path_cdfs_folder:
    :param n_threads:
    :return:
    """

    def exact_cdf(file_pw):
        """ Compute CDF for specific Multivatiate Bernoulli individuals (traget_indivs) """

        pw = np.loadtxt(file_pw)
        probs = pw[0]
        weights = pw[1]
        pmf, outcomes = exact_pmf(probs, weights)

        # Load target population of Multiv-Bern individuals
        # Compute their weighted outcomes (values) and corresponding CDFs
        file_pop = file_pw[:-len('pw.txt')] + 'pop.txt'
        target_indivs = np.loadtxt(file_pop)
        target_values = list(map(lambda x: np.dot(x, weights), target_indivs))
        target_cdfs = list(map(lambda value: sum(pmf[outcomes <= value]), target_values))

        file_cdfs = file_pw[:-len('pw.txt')] + 'exact.txt'

        np.savetxt(path_cdfs_folder + file_cdfs,
                   target_cdfs, fmt='%.16f')

    #path_disrt_folder = 'tests/data_test_precision/distr_d_weights_n0500/'

    files = glob.glob(path_disrt_folder + '*_pw.txt')
    with Pool(n_threads) as workers:
        pmap = workers.map
        pmap(exact_cdf, files)




if __name__ == '__main__':
    # generate distribution
    path_to_out_folder = 'tests/data_test_precision/tmp/'