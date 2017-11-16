import numpy as np
import bernmix_double.bernmix_double as bmd
import bernmix_double.bernmix_correction as bmc
from multiprocess import Pool
from itertools import count


def cdfs(i_distr):

    file_dir_pw = 'test_examples_n0100/'
    pw = np.loadtxt(file_dir_pw + 'distr' + str(i_distr).zfill(2) + '_pw.txt')
    probs = pw[0]
    weights = pw[1]

    target_indivs = np.loadtxt(file_dir_pw + 'distr' + str(i_distr).zfill(2) + '_pop.txt').astype(int)
    #cdf_corr = bmc.cdf_corrected(weights, probs, target_indiv, M, 1000)
    #cdfs_8 = bmd.cdf_fixed_scaling(probs, weights, target_indiv, 10 ** 4)
    #cdf_M = bmd.cdf_fixed_scaling(probs, weights, target_indiv, M)


    file_dir_res = 'res_examples_n0100/'

    for i in range(4, 9):
        if i == 4:
            n_repeats = 20
        elif i == 5:
            n_repeats = 15
        elif i == 6:
            n_repeats = 10
        elif i == 7:
            n_repeats = 3
        elif i == 8:
            n_repeats = 3
        print('Proc', i_distr, 'pow', i)
        cdfs = list(map(lambda x: np.array(list(bmd.cdfs_fixed_M(probs, weights, target_indivs, 10 ** i + x))),
               np.random.randint(0, 10 ** 3, n_repeats)))

        mx = np.empty((0, len(target_indivs)))
        for cdf in cdfs:
            mx = np.append(mx, [cdf], axis=0)

        np.savetxt(file_dir_res + 'distr' + str(i_distr).zfill(2) + '_M' + str(i) +'.txt', mx, fmt='%.16f')

    print('Proc', i_distr, 'is ended')

n_disrt = 20


with Pool(10) as workers:
    pmap = workers.map

    pmap(cdfs, range(n_disrt))








