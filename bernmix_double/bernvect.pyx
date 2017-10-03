cimport bernvect_evol

from libcpp.vector cimport vector

def evolving(pop, fitness, weights, mut_prob, n_epoch):
    # ANNA
    def copy_bern_vec(v):
        cdef vector[int] v_tmp
        for b in v:
            v_tmp.push_back(b)
        return v_tmp


    cdef int n_bern = len(weights)
    cdef vector[vector[int]] bernvectors
    cdef vector[double] cfitness
    for v in pop:
        bernvectors.push_back(copy_bern_vec(v))

    for f in fitness:
        cfitness.push_back(f)

    bernvect_evol.evolving(bernvectors, cfitness, weights, mut_prob, n_epoch)


    for i in range(bernvectors.size()):
        for j in range(bernvectors.at(i).size()):
            pop[i][j] = bernvectors.at(i).at(j)

    for i in range(cfitness.size()):
        fitness[i] = cfitness.at(i)

    return (pop, fitness)





print "Bernvector evolution is loaded"


