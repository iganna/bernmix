import numpy as np
from multiprocess import Pool



def cdfs(i_distr):

    def calc_prob(indiv, probs):
        prob_multiply = list(map(lambda i: probs[i] if indiv[i] == 1 else (1 - probs[i]),
                                 range(30)))
        return np.prod(prob_multiply)

    file_dir_pw = 'test_examples_n0030/'
    pw = np.loadtxt(file_dir_pw + 'distr' + str(i_distr).zfill(2) + '_pw.txt')
    probs = pw[0]
    weights = pw[1]

    target_indivs = np.loadtxt(file_dir_pw + 'distr' + str(i_distr).zfill(2) + '_pop.txt').astype(int)
    indiv_values = np.zeros(2 ** 30)
    indiv_probs = np.zeros(2 ** 30)

    indiv_values[0:2] = [0, weights[0]]
    indiv_probs[0] = calc_prob(np.zeros(30), probs)
    indiv_probs[1] = calc_prob(np.append([1], np.zeros(29)), probs)


    for i in range(1,30):
        print(i)
        n = 2 ** i
        indiv_values[n:2*n] = indiv_values[0:n] + weights[i]
        indiv_probs[n:2 * n] = indiv_probs[0:n] / (1 - probs[i]) * probs[i]



    print(sum(indiv_probs))
    # target_values = list(map(lambda x: np.dot(x, weights), target_indivs))
    #
    # cdfs = list(map(lambda x: sum(indiv_probs[indiv_values <= x]), target_values))
    #
    # file_dir_res = 'res_examples_n0030/'
    #
    # np.savetxt(file_dir_res + 'distr' + str(i_distr).zfill(2) + '_exact.txt',
    #            cdfs, fmt='%.16f')

    print('Proc', i_distr, 'is ended')

n_disrt = 20

with Pool(10) as workers:
    pmap = workers.map
    pmap(cdfs, range(n_disrt))


