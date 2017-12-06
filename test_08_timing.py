import numpy as np
import bernmix_int.bernmix_int as bmi
from scipy.stats import binom
import math
import time
import glob

from poibin import PoiBin
from multiprocess import Pool

# def MN(p1, p2):
#     return max(abs(p1-p2))
#
#
# def RSS(p1, p2):
#     return math.sqrt(sum((p1 - p2) ** 2))
#

#
#
# def accuracy_binomial():
#     for n in [10, 50, 100, 500, 1000, 5000]:
#
#         mn = 0
#         rss = 0
#         n_repeats = 20
#         for _ in range(n_repeats):
#             weights = np.ones(n)
#
#             # # Binomial distribution
#             # p = np.random.rand(1)
#             # probs = np.repeat(p, n)
#             # pmf_exact = binom.pmf(range(0, n + 1), n, p)
#
#             #PoissonBinomial distribution
#             probs = np.random.rand(n)
#             d = PoiBin(probs)
#             pmf_exact = d.pmf(range(0, n + 1))
#
#             pmf_bermnix = bmi.pmf(probs, weights)
#
#
#
#             mn += MN(pmf_bermnix, pmf_exact)
#             rss += RSS(pmf_bermnix, pmf_exact)
#             #print(mn, rss)
#         yield mn/n_repeats, rss/n_repeats
#         print(mn / n_repeats, rss / n_repeats)
#
#
# a = accuracy_binomial()
#
# list(a)


# #for n in [500,1000,5000,10000,25000,50000,100000,150000,200000,250000,300000,350000,400000,450000,500000,550000,650000,750000,900000,1000000]:
# for n in [500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000]:
#     weights = -np.ones(n)
#     n_repeat = 10
#     t1 = 0
#     t2 = 0
#     for _ in range(n_repeat):
#         probs = np.random.rand(n)
#
#         start1 = time.time()
#         d = PoiBin(probs)
#         pmf_exact = d.pmf(range(0, n + 1))
#
#         end1 = time.time()
#         t1 += end1 - start1
#
#
#         start2 = time.time()
#         pmf_bermnix = bmi.pmf(probs, weights)
#         end2 = time.time()
#         t2 += end2 - start2
#     t1 /= n_repeat
#     t2 /= n_repeat
#     print(t1, t2)


##======================================



def timing(file):

    pw = np.loadtxt(file)
    probs = pw[0]
    weights = pw[1]
    n = len(weights)

    start1 = time.time()
    d = PoiBin(probs)
    pmf_wieghted = bmi.pmf(probs, weights)
    end1 = time.time()
    t1 = end1 - start1

    weights = np.ones(n)
    start2 = time.time()
    pmf_unweighted = bmi.pmf(probs, weights)
    end2 = time.time()
    t2 = end2 - start2
    print(file[29:34], t1, t2)



files = glob.glob('timing/data_for_timing/*')
with Pool(10) as workers:
    pmap = workers.map
    pmap(timing, files)

