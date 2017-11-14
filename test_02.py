import bernmix_double.bernmix_double as bmd
import bernmix_double.gen_alg as ga
import numpy as np


n = 50
# weights = np.random.rand(n)
# probs = np.random.rand(n)
# target_indiv = np.random.randint(0, 2, n)
# np.savetxt("weights_n50.txt", weights, fmt='%.16f')
# np.savetxt("probs_n50.txt", probs, fmt='%.16f')
# np.savetxt("target_indiv_n50.txt", target_indiv, fmt='%i')


weights = np.loadtxt("weights_n50.txt", delimiter='\t')
probs = np.loadtxt("probs_n50.txt", delimiter='\t')
target_indiv = np.loadtxt("target_indiv_n50.txt", delimiter='\t').astype(int)

n_range = [10**6, 10**6+1000]

#
# cdf = []
# for m in range(n_range[0], n_range[1] + 1, 11):
#     cdf += [bmd.cdf_fixed_scaling(probs, weights, target_indiv, m)]

cdf_ethalon = 0.5653

m = 10**4

bmd.cdf_fixed_scaling(probs, weights, target_indiv, m)

k = 5
cdf, pmf = bmd.cdf_range(probs, weights, target_indiv, m, k)

mut_prob = 2 / n
pop_size = 1000
n_epoch = 500
target_value = np.dot(target_indiv, weights)
pop = ga.EvolvingPop(weights, pop_size, mut_prob)
pop.evolve(n_epoch, target_value)


probs_around_target_z = pop.single_m_overlap(m, target_indiv, probs, k)
s = sum(probs_around_target_z)
s[s == 0] = 1
probs_around_target_z[1] = probs_around_target_z[1] / s
probs_around_target_z[0] = probs_around_target_z[0] / s
probs_around_target_z


