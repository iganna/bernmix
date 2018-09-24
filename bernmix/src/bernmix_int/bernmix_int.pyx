cimport bernmix_fourier
from libc.stdlib cimport malloc, free

def pmf(probs, weights):
    # ANNA
    cdef int n = len(weights)
    cdef int sw = 0
    for wi in weights:
        sw += wi

    cdef int *w_new = <int *>malloc(n * sizeof(int))
    cdef double *p_new = <double *>malloc(n * sizeof(double))
    for i in range(n):
        w_new[i] = weights[i]
        p_new[i] = probs[i]

    cdef double *pmf_res = <double *> malloc((sw + 1) * sizeof(double))
    bernmix_fourier.pmf_bern_mixture(n, p_new, w_new, pmf_res)

    pmf_result = [ pmf_res[i] for i in range(sw+1) ]
    free(pmf_res)
    return pmf_result



