import numpy as np
import cvxopt
from cvxopt import glpk
import bernmix_double.bernmix_double as bmd
import bernmix_int.bernmix_int as bm_int


n = 100
w = np.random.rand(n)
probs = np.random.rand(n)
target_indiv = np.random.randint(0, 2, n).astype(int)



M = 10 ** 6


cdf_corr = bmc.cdf_corrected(w, probs, target_indiv, M)
cdf_6 = bmd.cdf_fixed_scaling(probs, w, target_indiv, 10 ** 6)
cdf_M = bmd.cdf_fixed_scaling(probs, w, target_indiv, M)
print(cdf_6, cdf_M, cdf_corr)



