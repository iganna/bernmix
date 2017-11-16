import numpy as np
import bernmix_int.bernmix_int as bm_int



def cdf_fixed_scaling(probs, weights, target_bern_vector, n_points):
    """
    This function calculates the CDF value for fixed number of points to approximate
    :param probs:
    :param weights:
    :param target_value:
    :param n_points:
    :return:
    """
    n = len(weights)
    w = sum(weights)
    c = (n_points + n) / w

    weights_new = np.around(weights * c, decimals=1).astype(int)
    new_target_value = np.dot(weights_new, target_bern_vector).astype(int)


    pmf = bm_int.pmf(probs, weights_new)


    return sum(pmf[:new_target_value + 1:])


def cdfs_fixed_M(probs, weights, target_vectors, n_points):
    """
    This function calculates the CDF value for fixed number of points to approximate
    :param probs:
    :param weights:
    :param target_value:
    :param n_points:
    :return:
    """
    n = len(weights)
    w = sum(weights)
    c = (n_points + n) / w

    weights_new = np.around(weights * c, decimals=1).astype(int)
    target_values = map(lambda x: np.dot(weights_new, x).astype(int), target_vectors)


    pmf = bm_int.pmf(probs, weights_new)

    pmfs = map(lambda x: sum(pmf[:x + 1:]), target_values)


    return pmfs



def cdf_range(probs, weights, target_bern_vector, n_points, k):
    """
    This function calculates the CDF value for fixed number of points to approximate
    :param probs:
    :param weights:
    :param target_value:
    :param n_points:
    :return:
    """
    n = len(weights)
    w = sum(weights)
    c = (n_points + n) / w

    weights_new = np.around(weights * c, decimals=1).astype(int)
    new_target_value = np.dot(weights_new, target_bern_vector)


    pmf = bm_int.pmf(probs, weights_new)


    return sum(pmf[:new_target_value + 1:]), \
           pmf[new_target_value - k: new_target_value + k + 1:]



def cdf_range_scaling(probs, weights, target_indiv, n_range,
                      n_epoch = 1000, pop_size = 1000, mut_prob = None):
    """
    This function calculates the CDF value for fixed number of points to approximate
    :param probs:
    :param weights:
    :param target_value:
    :param n_points:
    :return:
    """


    # initialisation of genetic algorithm
    mut_prob = 2/len(weights)
    pop = ga.EvolvingPop(weights, pop_size, mut_prob)

    target_value = pop.weighted_sum(target_indiv)

    pop.init_by_target(target_indiv)

    pop.evolve(n_epoch, target_value)

    # parameters for approximation

    m1 = n_range[0]
    m2 = n_range[1]
    n_points, indexes_down, indexes_up = pop.approximate(m1, m2, target_value)


    prob_discr = 0
    prob_discr = pop.prob_discrepancy(
                     n_points,
                     indexes_down,
                     indexes_up,
                     target_indiv,
                     probs)
    print(prob_discr)

    cdf = cdf_fixed_scaling(probs, weights, target_indiv, n_points)

    return  cdf, cdf + prob_discr

def cdf_range_scaling_not_corrected(probs, weights, target_indiv, n_range,
                      n_epoch = 1000, pop_size = 1000, mut_prob = None):
    """
    This function calculates the CDF value for fixed number of points to approximate
    :param probs:
    :param weights:
    :param target_value:
    :param n_points:
    :return:
    """


    # initialisation of genetic algorithm
    mut_prob = 2/len(weights)
    pop = ga.EvolvingPop(weights, pop_size, mut_prob)

    target_value = pop.weighted_sum(target_indiv)

    pop.evolve(n_epoch, target_value)

    # parameters for approximation

    m1 = n_range[0]
    m2 = n_range[1]
    n_points, indexes_down, indexes_up = pop.approximate(m1, m2, target_value)



    return cdf_fixed_scaling(probs, weights, target_indiv, n_points)



def cdf_permutations(probs, weights, target_indiv, n_permut):
    """
    Get CDF py permutations/simulations
    :param probs: 
    :param weights: 
    :param target_indiv: 
    :param n_permut: 
    :return: 
    """

    pop = np.empty(shape=[0, n_permut], dtype=int)
    for i in range(len(weights)):
        pop = np.append(pop, np.random.binomial(1, probs[i], (1, n_permut)), axis=0)

    pop = np.transpose(pop)

    target_value = np.dot(weights, target_indiv)

    pop_values = list(map(lambda x: np.dot(x, weights), pop))

    cdf_target_value = sum(value <= target_value for value in pop_values)/n_permut

    return cdf_target_value







