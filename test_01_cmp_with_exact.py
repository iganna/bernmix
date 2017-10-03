import bernmix_double.bernmix_double as bmd
import bernmix_double.gen_alg as ga
import numpy as np


n = 20
weights = np.random.rand(n)
probs = np.random.rand(n)
pop_size = 1000
mut_prob = 0.2
n_epoch = 100


target_indiv = np.random.randint(0, 2, n)

M = 10**5
n_range = (M - 100, M)
bmd.cdf_range_scaling(probs, weights, target_indiv, n_range,
                      n_epoch, pop_size, 2/len(weights))


n_permut = 10000000
bmd.cdf_permutations(probs, weights, target_indiv, n_permut)
#----------------------------------------------
# initialisation
weights = np.random.rand(20)
pop_size = 1000
mut_prob = 0.2
n_epoch = 300

pop = ga.EvolvingPop(weights, pop_size, mut_prob)

pop_init = pop.pop
target_indiv = pop.rand_indiv()
target_value = pop.weighted_sum(target_indiv)
# target_value = np.random.random() * sum(weights)

pop.evolve(n_epoch, target_value)
sum(pop.fitness)

pop.pop = pop_init
pop.evolve_new(n_epoch, target_value)
sum(pop.fitness)


# parameters for approximation
m2 = 100000
m1 = m2 - 1000
c, x, possible_indiv = pop.approximate(m1, m2, target_value)

