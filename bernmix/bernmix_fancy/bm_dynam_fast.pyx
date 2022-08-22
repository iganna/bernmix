# print("bernmix dynamic alg")
from libc.stdlib cimport malloc, free

def dp_int_fast(p, w):
    n = len(w)
    cdef w_tot = sum(w)
    cdef w_curr = 1

    cdef double *s = <double *>malloc((w_tot+1) * sizeof(double))
    cdef double *s_tmp = <double *>malloc((w_tot+1) * sizeof(double))

    s[0] = 1

    for i in range(n):

        for j in range(w[i]+1):
            s_tmp[j] = s[j] * (1 - p[i])
        for j in range(w[i], w_curr):
            s_tmp[j] = s[j] * (1 - p[i]) + s[j- w[i]] * p[i]
        for j in range(w_curr, w_curr + w[i]):
            if(j- w[i] >= 0):
                s_tmp[j] = s[j- w[i]] * p[i]
            else:
                s_tmp[j] = 0

        s, s_tmp = s_tmp, s
        w_curr += w[i]

        #s_print = [s[k] for k in range(w_curr)]
        #print(s_print)


    s_result = [ s[i] for i in range(w_tot+1) ]
    free(s_tmp)
    free(s)

    return s_result
