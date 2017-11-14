import numpy as np

import bernmix_double.bernmix_double as bmd
import bernmix_double.bernmix_correction as bmc
import bernmix_int.bernmix_int as bm_int


def prob_of_indiv(indiv, prob):
    n = len(indiv)
    prob_multiply = list(map(lambda i: prob[i] if indiv[i] == 1 else (1 - prob[i]),
                             range(n)))
    return np.prod(prob_multiply)

# initialise parameters weights and probs

file_name_pw = '/Users/anna/OneDrive/polytech/stud_projects/enrichment/scripts/bernmix/test_examples_n=27/ex1_n27_pw'
mx_pw = np.loadtxt(file_name_pw, skiprows=1, delimiter='\t')
file_name_points = '/Users/anna/OneDrive/polytech/stud_projects/enrichment/scripts/bernmix/test_examples_n=27/ex1_n27_points'
mx_points = np.loadtxt(file_name_points, skiprows=1, delimiter='\t')
probs = mx_pw[:, 0]
weights = mx_pw[:, 1]
n = len(weights)




n_repeats = 5
cdf_perm = np.zeros((n_repeats, 9))
for i_perm in range(n_repeats):
    print(i_perm)
    # get random point
    point_info = mx_points[np.random.randint(0, len(mx_points), 1),:]
    target_indiv = point_info[0, 0:n]
    target_value = np.dot(target_indiv, weights)
    target_cdf = point_info[0, n+1]


    # by permutations
    n_permutations = 1000000
    cdf_perm[i_perm, 0] = abs(bmd.cdf_permutations(probs, weights, target_indiv, n_permutations) - target_cdf)

    idx_cdf_next = 1
    # by our method
    for M in map(lambda x: 10 ** x, range(4, 8)):
        print(M)
        cdf_perm[i_perm, idx_cdf_next] = abs(bmd.cdf_fixed_scaling(probs, weights, target_indiv, M) - target_cdf)
        idx_cdf_next += 1
        cdf_perm[i_perm, idx_cdf_next] = abs(bmc.cdf_corrected(weights, probs, target_indiv, M) - target_cdf)
        idx_cdf_next += 1

