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


    file_dir_res = 'res_examples_n0100/'
    n_solutions = 100

    for i in range(4, 7):

        print('Proc', i_distr, 'pow', i)

        cdfs = list(map(lambda target_indiv: bmc.cdf_corrected(weights, probs, target_indiv, 10 ** i, n_solutions), target_indivs))

        np.savetxt(file_dir_res + 'distr' + str(i_distr).zfill(2) + '_M' + str(i) +
                   '_corrected' + str(n_solutions).zfill(4) + '.txt',
                   [np.array(cdfs)], fmt='%.16f')

    print('Proc', i_distr, 'is ended')

n_disrt = 20


with Pool(10) as workers:
    pmap = workers.map
    pmap(cdfs, range(n_disrt))








