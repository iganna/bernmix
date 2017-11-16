
import bernmix_double.bernmix_double as bmd
import bernmix_double.gen_alg as ga
import numpy as np

# ===============================
# generate distributions
# 5 distributions, 20 points in each
# ===============================
n = 100
pop_size = 5
n_disrt = 20

for i_distr in range(n_disrt):
    w = np.random.rand(n)
    p = np.random.rand(n)
    pw = np.concatenate(([p], [w]), axis=0)
    file_dir_pw = '/Users/anna/OneDrive/polytech/stud_projects/enrichment/scripts/bernmix/test_examples_n0100/'
    np.savetxt(file_dir_pw + 'distr' + str(i_distr).zfill(2) + '_pw.txt', pw, fmt='%.16f')


    pop = np.random.randint(0, 2, (pop_size, n))
    np.savetxt(file_dir_pw + 'distr' + str(i_distr).zfill(2) + '_pop.txt', pop, fmt='%i')





# ===============================
# generate populations for 100 random points:
# 5 distributions, 20 points in each
# ===============================

n = 100
pop_size = 5
n_epoch = 500
M = 10**4
n_range = (M - 1000, M)

i_distr = 0
i_point = 0

for i_distr in range(5):
    file_name_pw = '/Users/anna/OneDrive/polytech/stud_projects/enrichment/scripts/bernmix/test_examples2/ex' + str(i_distr+1) + '_n20_pw'
    mx_pw = np.loadtxt(file_name_pw, skiprows=1, delimiter='\t')
    file_name_points = '/Users/anna/OneDrive/polytech/stud_projects/enrichment/scripts/bernmix/test_examples2/ex' + str(i_distr+1) + '_n20_points'
    mx_points = np.loadtxt(file_name_points, skiprows=1, delimiter='\t')
    probs = mx_pw[:, 0]
    weights = mx_pw[:, 1]
    for i_point in range(25):
        print(i_distr+1, i_point+1)
        target_indiv = mx_points[i_point, 0:n].astype(int)
        target_value = np.dot(target_indiv, weights)
        print(target_value - mx_points[i_point,n])
        mut_prob = 4 / n
        pop = ga.EvolvingPop(weights, pop_size, mut_prob)
        pop.init_by_target(target_indiv)
        pop.evolve(n_epoch, target_value)

        for M in range(10**4, 10**4+201,11):
            print(pop.prop_discrepancy(M+2, target_indiv, probs))


        file_population = '/Users/anna/OneDrive/polytech/stud_projects/enrichment/scripts/bernmix/test_output_pop_0500/ex' + str(i_distr+1) + '_' + str(i_point+1) + '.txt'
        np.savetxt("pop_tmp", pop.pop, fmt='%i')

# ===============================
# Make the figure: (m, n_overlap)
# ===============================


weights = np.random.randint(0,200,20)

i_distr = 0
i_point = 0

for i_distr in range(5):
    file_name_pw = '/Users/anna/OneDrive/polytech/stud_projects/enrichment/scripts/bernmix/test_examples2/ex' + str(i_distr+1) + '_n20_pw'
    mx_pw = np.loadtxt(file_name_pw, skiprows=1, delimiter='\t')
    file_name_points = '/Users/anna/OneDrive/polytech/stud_projects/enrichment/scripts/bernmix/test_examples2/ex' + str(i_distr+1) + '_n20_points'
    mx_points = np.loadtxt(file_name_points, skiprows=1, delimiter='\t')
    probs = mx_pw[:, 0]
    weights = mx_pw[:, 1]
    for i_point in range(25):
        file_population = '/Users/anna/OneDrive/polytech/stud_projects/enrichment/scripts/bernmix/test_output_pop_0500/ex' + str(i_distr + 1) + '_' + str(i_point + 1) + '.txt'
        print(i_distr+1, i_point+1)
        target_indiv = mx_points[i_point, 0:n].astype(int)
        target_value = np.dot(target_indiv, weights)
        print(target_value - mx_points[i_point,n])
        mut_prob = 2 / n
        pop = ga.EvolvingPop(weights, pop_size, mut_prob)

        pop.evolve(n_epoch, target_value)

        # pop.pop = np.loadtxt(file_population)

        n_overlaps, d_overlaps, rnage_onerlaps, descr = pop.M_overlap(n_range[0], n_range[1], target_indiv, probs)




        res = np.transpose(np.concatenate((np.matrix(n_overlaps),
                              np.matrix(d_overlaps),
                              np.matrix(rnage_onerlaps),
                              np.matrix(descr)
                              ), axis=0))


        np.savetxt("res_uniform.txt", res, fmt='%.16f')

        # mx_points[i_point, n+1]
        #
        cdf = []
        for m in range(n_range[0], n_range[1]+1):
            cdf += [bmd.cdf_fixed_scaling(probs, weights, target_indiv, m)]


        np.savetxt("cdf_uniform.txt", cdf, fmt='%.16f')


        s = "res_" + str(i_distr) + '_' + str(i_point)
        np.savetxt(s, res, fmt='%.16f')

#===================================================
# tmp
#===================================================





n = 50
weights = np.random.rand(n)
probs = np.random.rand(n)
n_range = [10**6, 10**6+1000]
target_indiv = np.random.randint(0, 2, n)

cdf = []
for m in range(n_range[0], n_range[1] + 1, 11):
    cdf += [bmd.cdf_fixed_scaling(probs, weights, target_indiv, m)]

np.savetxt("cdf_n100.txt", cdf, fmt='%.16f')




target_indiv = np.random.randint(0, 2, n)

# read weights and probs from file
file_name_pw = '/Users/anna/OneDrive/polytech/stud_projects/enrichment/scripts/bernmix/test_examples2/ex1_n20_pw'
mx_pw = np.loadtxt(file_name_pw, skiprows=1, delimiter='\t')
file_name_points = '/Users/anna/OneDrive/polytech/stud_projects/enrichment/scripts/bernmix/test_examples2/ex1_n20_points'
mx_points = np.loadtxt(file_name_points, skiprows=1, delimiter='\t')


probs = mx_pw[:,0]
weights = mx_pw[:,1]
target_indiv = mx_points[2,0:n].astype(int)

target_value = np.dot(target_indiv, weights)
target_value
mx_points[2,n]


bmd.cdf_range_scaling(probs, weights, target_indiv, n_range,
                      n_epoch, pop_size, 2/len(weights))

bmd.cdf_fixed_scaling(probs, weights, target_indiv, M)


mx_points[2,n+1]

n_permut = 100000
bmd.cdf_permutations(probs, weights, target_indiv, n_permut)
