import numpy as np
import glob
from multiprocess import Pool



def pmfs(i_distr):

    def calc_prob(indiv, probs, n):
        prob_multiply = list(map(lambda i: probs[i] if indiv[i] == 1 else (1 - probs[i]),
                                 range(n)))
        return np.prod(prob_multiply)


    file_dir_pw = 'test_int_weights/'

    for file in glob.glob(file_dir_pw + '*' + str(i_distr).zfill(2) + '_pw.txt'):

        pw = np.loadtxt(file)
        probs = pw[0]
        weights = np.array(pw[1]).astype(int)
        n = len(probs)
        W = sum(weights)


        indiv_values = np.zeros(2 ** n)
        indiv_probs = np.zeros(2 ** n)
        indiv_values[0:2] = [0, weights[0]]
        indiv_probs[0] = calc_prob(np.zeros(n), probs, n)
        indiv_probs[1] = calc_prob(np.append([1], np.zeros(n-1)), probs, n)
        for i in range(1,n):
            print(i)
            m = 2 ** i
            indiv_values[m:2 * m] = indiv_values[0:m] + weights[i]
            indiv_probs[m:2 * m] = indiv_probs[0:m] / (1 - probs[i]) * probs[i]

        print('Proc', file[17:32], 'pmf is computed')

        pmf_bm = np.zeros(W + 1)
        for i in range(0, W+1):
            pmf_bm[i] = sum(indiv_probs[indiv_values == i])


        file_dir_res = 'res_int_weights/'

        np.savetxt(file_dir_res +
                   file[17:32] + '_exact_pmf.txt',
                   pmf_bm, fmt='%.16f')

        print('Proc', file[17:32], 'is ended')

n_disrt = 20

with Pool(10) as workers:
    pmap = workers.map
    pmap(pmfs, range(n_disrt))



