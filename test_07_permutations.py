import numpy as np
import bernmix_double.bernmix_double as bmd
from multiprocess import Pool


def cdfs_calc(i_distr):

    file_dir_pw = 'test_examples_n0100/'
    pw = np.loadtxt(file_dir_pw + 'distr' + str(i_distr).zfill(2) + '_pw.txt')
    probs = pw[0]
    weights = pw[1]

    target_indivs = np.loadtxt(file_dir_pw + 'distr' + str(i_distr).zfill(2) + '_pop.txt').astype(int)


    file_dir_res = 'res_examples_n0100/'
    n_permut = 1000000


    cdfs = list(map(lambda target_indiv: bmd.cdf_permutations(probs, weights, target_indiv, n_permut), target_indivs))

    np.savetxt(file_dir_res + 'distr' + str(i_distr).zfill(2)  +
               '_permut' + str(n_permut).zfill(4) + '.txt',
               [np.array(cdfs)], fmt='%.16f')



n_disrt = 20


with Pool(10) as workers:
    pmap = workers.map
    pmap(cdfs_calc, range(n_disrt))

