cdef extern from "bernmix_fourier.h":
    void pmf_bern_mixture(int n, double* P, int* W, double* res)
