from libcpp.vector cimport vector
cdef extern from "bernvect_evol.h":
    void evolving (vector[vector[int]] & bernvectors,
                   vector[double] & fitness,
                   vector[double] weights,
                   double mut_prob,
                   int n_epoch);
