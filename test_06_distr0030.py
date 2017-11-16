import numpy as np





def cdfs(i_distr):

    def val_prob(indiv, probs, weights):
        prob_multiply = list(map(lambda i: probs[i] if indiv[i] == '1' else (1 - probs[i]),
                                 range(30)))
        val_multiply = list(map(lambda i: weights[i] if indiv[i] == '1' else 0,
                                 range(30)))
        return np.prod(prob_multiply), sum(val_multiply)


    file_dir_pw = 'test_examples_n0030/'
    pw = np.loadtxt(file_dir_pw + 'distr' + str(i_distr).zfill(2) + '_pw.txt')
    probs = pw[0]
    weights = pw[1]

    target_indivs = np.loadtxt(file_dir_pw + 'distr' + str(i_distr).zfill(2) + '_pop.txt').astype(int)

    indiv_values = np.zeros(2 ** 30)
    indiv_probs = np.zeros(2 ** 30)
    for indiv_id in range(0, 2 ** 30):
        indiv = bin(indiv_id)[2:]
        indiv += '0' * (30 - len(indiv))
        indiv_probs[indiv_id], indiv_values[indiv_id] = val_prob(indiv, probs, weights)

    target_values = list(map(lambda x: np.dot(x, weights), target_indivs))

    cdfs = list(map(lambda x: sum(indiv_probs[indiv_values <= x]), target_values))

    file_dir_res = 'res_examples_n0030/'

    np.savetxt(file_dir_res + 'distr' + str(i_distr).zfill(2) + '_exact.txt',
               cdfs, fmt='%.16f')

    print('Proc', i_distr, 'is ended')

n_disrt = 20

with Pool(10) as workers:
    pmap = workers.map
    pmap(cdfs, range(n_disrt))



