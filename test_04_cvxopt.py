import numpy as np
import bernmix_double.bernmix_double as bmd
import bernmix_double.bernmix_correction as bmc
from multiprocess import Pool
from itertools import count


def cdfs(a):
    run_id, w, probs = a
    target_indiv = np.random.randint(0, 2, n).astype(int)
    cdf_corr = bmc.cdf_corrected(w, probs, target_indiv, M)
    cdf_6 = bmd.cdf_fixed_scaling(probs, w, target_indiv, 10 ** 6)
    cdf_M = bmd.cdf_fixed_scaling(probs, w, target_indiv, M)
    print('Proc ', run_id, ' is ended')
    return cdf_6, abs(cdf_M - cdf_6), abs(cdf_corr - cdf_6)


n = 100
M = 10 ** 4


#with Pool(3) as workers:
    # pmap = workers.map

p = Pool(3)
pmap = p.map
x = pmap(cdfs, zip(count(), [np.random.rand(n) for _ in range(10)], [np.random.rand(n) for _ in range(10)]))


for xi in x:
    print(xi)





