"""

This module generates data for verifying run-time complexity of PMF calculation
and measures the run-time

"""

__author__ = "Anna Igolkina"
__all__ = ["gen_int_distr", "time_measuring"]

import time
import glob
import os
from multiprocess import Pool
import numpy as np
import itertools as it
import bernmix_int.bernmix_int as bmi


def ten_in_power(powers):
    return [10 ** x for x in powers]


def gen_int_distr(grid_for_n, grid_for_weight_sum,
                               n_repeats, path_to_out_folder):
    """
    This function generates parameters of distributions of positively-weighted sum of BRVs
    with integer weights, which sum equals to the weight_sum
    :param grid_for_n:
    :param grid_for_weight_sum:
    :param n_repeats:
    :param path_to_out_folder:
    :return:
    """
    cnt = 0
    n_digits = len(str(len(grid_for_n) * len(grid_for_weight_sum) * n_repeats))
    for n, weight_sum, _ in it.product(grid_for_n, grid_for_weight_sum, range(n_repeats)):

        # cannot find n positive integer numbers which sum to weight_sum
        if weight_sum < n:
            continue

        # randomly define probabilities of i.ni.d. BRVs
        probs = np.random.rand(n)

        # define weights
        weights = np.random.rand(n)
        weights = np.round(weights * weight_sum / sum(weights)).astype(int)

        # no weights should be equal to zero
        weights[weights == 0] += 1

        # remove discrepansy between the sum of weights and the target sum
        while weight_sum != sum(weights):
            if weight_sum > sum(weights):
                weights[np.random.randint(n)] += 1
            else:
                weights[weights == max(weights)] -= 1

        # save probabilities and weights into the file
        np.savetxt(path_to_out_folder + str(cnt).zfill(n_digits) + '.txt',
                   np.concatenate(([probs], [weights]), axis=0), fmt='%0.6f')
        cnt += 1

        yield '{}\t{}\t{}'.format(cnt, n, weight_sum)


def time_measuring(path_to_folder_distr, path_to_folder_time, n_threads):
    """

    :param path_to_folder_distr:
    :param path_to_folder_time:
    :param n_threads: The number of worker threads to execute
    :return:
    """

    def timing(file):

        # read
        pw = np.loadtxt(file)
        probs = pw[0]
        weights = pw[1]

        # time measuring
        start1 = time.time()
        bmi.pmf(probs, weights)
        end1 = time.time()
        t1 = end1 - start1

        # save
        distr_name = file.split('/')[2]
        np.savetxt(path_to_folder_time + 'time_' + distr_name, [t1], fmt='%0.8f')

    # create out folder if it does not exist
    if not os.path.exists(path_to_folder_time):
        os.makedirs(path_to_folder_time)

    files = glob.glob(path_to_folder_distr + 'distr*.txt')
    with Pool(n_threads) as workers:
        pmap = workers.map
        pmap(timing, files)


if __name__ == '__main__':

    path_to_folder_distr = 'tests/data_test_runtime/'
    path_to_folder_time = 'tests/res_test_runtime/'
    #
    # n_repeats = 10
    #
    # # generate distributions when n varies while sum of weights is fixed
    # grid_for_n = np.round(ten_in_power(np.arange(1, 4.1, 0.2))).astype(int)
    # grid_for_weight_sum = np.round(ten_in_power(np.arange(3, 6.1, 1))).astype(int)
    # path_to_out_folder = path_to_folder_distr + 'distr_var_n_'
    # info = gen_int_distr(grid_for_n, grid_for_weight_sum,
    #                      n_repeats, path_to_out_folder)
    # np.savetxt(path_to_folder_distr + 'info_var_n.csv',
    #            list(info), '%s')
    #
    # # generate distributions when n is fixed while sum of weights varies
    # grid_for_n = np.round(ten_in_power(np.arange(2, 4.1, 1))).astype(int)
    # grid_for_weight_sum = np.round(ten_in_power(np.arange(3, 6.1, 0.2))).astype(int)
    # path_to_out_folder = path_to_folder_distr + 'distr_var_w_'
    # info = gen_int_distr(grid_for_n, grid_for_weight_sum,
    #                      n_repeats, path_to_out_folder)
    # np.savetxt(path_to_folder_distr + 'info_var_w.csv',
    #            list(info), '%s')

    # timing
    n_threads = 10
    time_measuring(path_to_folder_distr, path_to_folder_time, n_threads)

