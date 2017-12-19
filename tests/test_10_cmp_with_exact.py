import numpy as np
import glob
import bernmix_int.bernmix_int as bmi
import math


def MN(p1, p2):
    return max(abs(p1-p2))


def RSS(p1, p2):
    return math.sqrt(sum((p1 - p2) ** 2))



file_dir_pw = 'test_int_weights/'
file_dir_res = 'res_int_weights/'

w = 5
x = np.zeros(6)
id = 0
for file_pw in glob.glob(file_dir_pw + '*' + '_pw.txt'):
    distr_name = file_pw[17:32]
    if int(distr_name[11]) != w:
        continue

    id += 1
    file_res = file_dir_res + distr_name + '_exact_pmf.txt'

    pw = np.loadtxt(file_pw)
    probs = pw[0]
    weights = np.array(pw[1]).astype(int)
    n = len(probs)

    pmf_exact = np.loadtxt(file_res)

    pmf_bm = bmi.pmf(probs, weights)

    mn = MN(pmf_exact, pmf_bm)
    rss = RSS(pmf_exact, pmf_bm)

    if n == 10:
        x[0] += mn
        x[1] += rss
    elif n == 20:
        x[2] += mn
        x[3] += rss
    else:
        x[4] += mn
        x[5] += rss

x /= 20
x






