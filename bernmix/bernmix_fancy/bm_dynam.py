import numbers
import numpy as np
import scipy.stats as stats
from scipy.stats import mvn
from scipy.optimize import minimize
from functools import partial
from random import seed
from itertools import product, combinations

def dp_int(p, w):
    """
    The dynamic algorithm to estimate PMF (probability mass function)
    :param p: vector of probabilities
    :param w: vector of weights
    :return: PMF
    """
    # Checks
    if len(p) != len(w):
        raise ValueError('Probabilities and weights should be of the same length')
    n = len(p)
    # print(n)
    if n == 0:
        raise ValueError('Probabilities and weights are empty')
    if not all([isinstance(w_tmp, int) for w_tmp in w]):
        raise ValueError('Weights are not integer')
    if not all([isinstance(p_tmp, numbers.Number) for p_tmp in p]):
        raise ValueError('Probabilities are not numbers')
    if not all([0 <= p_tmp <= 1 for p_tmp in p]):
        raise ValueError('Probabilities are not in [0, 1]')

    # Algorithm
    w_tot = sum(w)
    w_curr = 1
    s = [1] + [0] * w_tot  # you don't need + 1 here
    for i in range(n):
        w_curr = w_curr + w[i]
        # s_tmp = [0] * w[i] + s

        # print([s_tmp[2], s[2 - w[i]]])
        s[0:(w_curr)] = [s[j] * (1 - p[i]) + s[j- w[i]] * p[i]
                         if j >= w[i]
                         else s[j] * (1 - p[i])
                            for j in range(w_curr)]
    return s



def prob_stat(pa):
    """
    For a given probability distribution, this function calculates
    Mean, Variabce and all quantiles (including -Inf and +Inf)
    :param pa: PMF
    :return:
    """
    ma = sum([i * p for i, p in enumerate(pa)])
    va = sum([(i - ma) ** 2 * p for i, p in enumerate(pa)])
    pc_cum = [1 if p > 1 else p for p in np.cumsum(pa)]

    points = [-np.inf] + list(stats.norm.ppf(pc_cum, loc=0, scale=1))
    return ma, va, points


def discrepancy(pa, pb, v_ab, cov_ab):
    """
    Returns the discrepancy between true covariance value (v_ab)
     and covariance value from bvNormal approximation with (cov_ab)
    :param pa: PMF of one variable
    :param pb: PMF of another variable
    :param v_ab: True value of covariance
    :param cov_ab: Tested value of covariance
    :return: discrepancy
    """

    if abs(cov_ab) > 1: return np.inf
    ma, va, points_a = prob_stat(pa)

    mb, vb, points_b = prob_stat(pb)

    cov_mx = [[1, cov_ab], [cov_ab, 1]]

    joint_distr = np.empty((len(pa), len(pb)))
    for i in range(len(pa)):
        for j in range(len(pb)):
            res = mvn.mvnun([points_a[i], points_b[j]],
                            [points_a[i + 1], points_b[j + 1]], [0, 0], cov_mx)

            joint_distr[i, j] = res[0]


    cov_new = 0
    for i in range(len(pa)):
        for j in range(len(pb)):
            cov_new += joint_distr[i, j] * (i - ma) * (j - mb)

    return abs(cov_new - v_ab)



def joint_distr_cov(pa, pb, v_ab, echo=False):
    """
    Estimate table of joint distribution between two RVs
    :param pa: PMF of the first RV
    :param pb: PMF of the second RV
    :param v_ab: required covariance between RVs
    :return: table with joint distribution
    """

    if abs(v_ab) >= 1e-5:
        f_opt = partial(discrepancy, pa, pb, v_ab)
        # bnds = ((-1, 1), ())
        seed(239)
        res = minimize(f_opt, 0, method='SLSQP', bounds=[(-1, 1)],
                       options={'ftol': 1e-8, 'disp': False})
        cov_ab = res.x
        if echo:
            print(f'Optimised covariance is {cov_ab}')
        # print(f'Optimised covariance is {cov_ab}')

        if(np.isnan(cov_ab)):
            if v_ab > 0:
                cov_ab = 1
            else:
                cov_ab = -1
            if echo:
                print(f'New is {cov_ab}')
    else:
        cov_ab = 0

    ma, va, points_a = prob_stat(pa)
    mb, vb, points_b = prob_stat(pb)

    cov_mx = [[1, cov_ab], [cov_ab, 1]]


    joint_distr = np.empty((len(pa), len(pb)))
    for i in range(len(pa)):
        for j in range(len(pb)):
            res = mvn.mvnun([points_a[i], points_b[j]],
                            [points_a[i + 1], points_b[j + 1]], [0, 0], cov_mx)

            joint_distr[i, j] = res[0]
            if(np.isnan(joint_distr[i, j])):
                joint_distr[i, j] = 0
            if echo:
                print([i, j, res])

    cov_new = 0
    for i in range(len(pa)):
        for j in range(len(pb)):
            cov_new += joint_distr[i, j] * (i - ma) * (j - mb)

    if echo:
        print([v_ab, cov_new])
        print(f'Discrepancy: {sum(sum(joint_distr) - pb)}, {sum(joint_distr[:, 1] + joint_distr[:, 0] - pa)}, {v_ab- cov_new}')

    return joint_distr


def est_distr(p, w, v, v_zero=False):
    """
    Estimate PMF for a wighted sum of BRVs with known covariance matrix between RVs
    :param p: probabilities
    :param w: weights
    :param v: covariance matrix between BRVs
    :param v_zero: To set all covariance az zero
    :return:
    """
    # ---------------------------------------------------
    # Checks
    # ---------------------------------------------------
    if len(p) != len(w):
        raise ValueError('Probabilities and weights should be of the same length')
    n = len(p)
    if n == 0:
        raise ValueError('Probabilities and weights are empty')
    if not all([isinstance(w_tmp, int) for w_tmp in w]):
        raise ValueError('Weights are not integer')
    if not all([isinstance(p_tmp, numbers.Number) for p_tmp in p]):
        raise ValueError('Probabilities are not numbers')
    if not all([0 <= p_tmp <= 1 for p_tmp in p]):
        raise ValueError('Probabilities are not in [0, 1]')

    # ---------------------------------------------------
    #       Algorithm
    # ---------------------------------------------------

    w_tot = sum(w)
    w_curr = 1
    s = [1] + [0] * w_tot  # you don't need + 1 here
    for i in range(n):
        if i == 0:
            w_curr = w_curr + w[i]
            s_tmp = [0] * w[i] + s
            s[0:(w_curr)] = [s[j] * (1 - p[i]) + s_tmp[j] * p[i] for j in range(w_curr)]
        else:
            # break
            pa = s[:w_curr]
            # print(pa)
            v_ab = sum([w[j] * v[i, j] for j in range(i)])
            # print(v_ab)
            if v_zero:
                v_ab = 0

            if np.isnan(v_ab):
                raise ValueError('Variance is nan')
            pb = [1 - p[i], p[i]]

            # print(pa)
            # print(pb)
            joint_distr = joint_distr_cov(pa, pb, v_ab)
            # print(joint_distr)

            s_tmp1 = [0] * w[i] + list(joint_distr[:, 1])
            s_tmp2 = list(joint_distr[:, 0]) + [0] * w[i]
            w_curr = w_curr + w[i]
            s[0:(w_curr)] = [x1 + x2 for x1, x2 in zip(s_tmp1, s_tmp2)]

            # print(s[0:(w_curr)])
    return(s)


def get_v_max(pa, pb):
    """
    Function to generate max value of variance
    :param pa:
    :param pb:
    :return:
    """

    f_opt = lambda cov_ab: (-1) * discrepancy(pa, pb, 0, cov_ab)

    res = minimize(f_opt, 0, method='SLSQP',
                   options={'ftol': 1e-8, 'disp': False})
    v_ab = (-1) * res.fun

    j_distr = joint_distr_cov(pa, pb, v_ab, echo=False)

    return v_ab, j_distr


def dyn_tree(p, w):
    if len(p) != len(w):
        raise ValueError('Probabilities and weights should be of the same length')
    n = len(p)
    print(n)
    if n == 0:
        raise ValueError('Probabilities and weights are empty')
    if not all([isinstance(w_tmp, int) for w_tmp in w]):
        raise ValueError('Weights are not integer')
    if not all([isinstance(p_tmp, numbers.Number) for p_tmp in p]):
        raise ValueError('Probabilities are not numbers')
    if not all([0 <= p_tmp <= 1 for p_tmp in p]):
        raise ValueError('Probabilities are not in [0, 1]')

    s_blocks = []
    for i in range(n):
        s_tmp = [1 - p[i]] + [0] * (w[i] - 1) + [p[i]]
        s_blocks += [s_tmp]

    while len(s_blocks) > 1:
        s_len = len(s_blocks)
        m = int((s_len - s_len % 2) / 2)
        s_blocks_new = []
        for j in range(m):
            i = j * 2
            # if i == s_len:
            #     s_blocks_new += [s_blocks[i]]
            #     break
            b1 = s_blocks[i]
            b2 = s_blocks[i + 1]
            b_new = [0] * (len(b1) + len(b2) - 1)
            for i1 in range(len(b1)):
                if b1[i1] == 0:
                    continue
                for i2 in range(len(b2)):
                    i_new = i1 + i2
                    b_new[i_new] += b1[i1] * b2[i2]
            s_blocks_new += [b_new]
        if s_len % 2 == 1:
            s_blocks_new += [s_blocks[s_len - 1]]
        s_blocks = s_blocks_new
    return s_blocks



def dyn_tree_cov(p, w, v):
    if len(p) != len(w):
        raise ValueError('Probabilities and weights should be of the same length')
    n = len(p)
    if n == 0:
        raise ValueError('Probabilities and weights are empty')
    if not all([isinstance(w_tmp, int) for w_tmp in w]):
        raise ValueError('Weights are not integer')
    if not all([isinstance(p_tmp, numbers.Number) for p_tmp in p]):
        raise ValueError('Probabilities are not numbers')
    if not all([0 <= p_tmp <= 1 for p_tmp in p]):
        raise ValueError('Probabilities are not in [0, 1]')

    s_blocks = []
    for i in range(n):
        s_tmp = [1 - p[i]] + [0] * (w[i] - 1) + [p[i]]
        s_blocks += [s_tmp]

    v_blocks = v.copy()
    for i, j in combinations(range(n), 2):
        v_blocks[i, j] *= w[i] * w[j]
        v_blocks[j, i] *= w[i] * w[j]

    while len(s_blocks) > 1:
        n = len(s_blocks)

        var_block = [i for _, i, _ in [prob_stat(s) for s in s_blocks]]

        v_norm = v_blocks.copy()
        for i, j in combinations(range(n), 2):
            v_norm[i, j] /= (var_block[i] * var_block[j]) ** (1 / 2)
            v_norm[j, i] /= (var_block[i] * var_block[j]) ** (1 / 2)


        # find the pair with the highest variance
        i = j = 0
        if v_norm.max() == 0:
            i = 0
            j = 1
        else:
            while i == j:
                v_norm[i, j] = 0
                k = np.argmax(v_norm)
                j = k % n
                i = round((k - j) / n)
        # print(i, j)
        # print([v_blocks[i, j], v_norm[i, j]])


        pa = s_blocks[i]
        pb = s_blocks[j]
        v_ab = v_blocks[i, j]

        # print(pa)
        # print(pb)
        j_distr = joint_distr_cov(pa, pb, v_ab)

        v_bb = 0

        b_new = [0] * (len(pa) + len(pb) - 1)
        for i1, i2 in product(range(len(pa)),
                                   range(len(pb))):
            b_new[i1 + i2] += j_distr[i1][i2]

        if i > j:
            s_blocks.pop(i)
            s_blocks.pop(j)
        else:
            s_blocks.pop(j)
            s_blocks.pop(i)
        s_blocks += [b_new]

        v_new = list(v_blocks[i] + v_blocks[j])
        v_blocks = np.append(v_blocks, [v_new], axis=0)
        v_blocks = np.append(v_blocks.transpose(), [v_new+ [0]], axis=0)
        v_blocks = np.delete(v_blocks, [i, j], axis=0)
        v_blocks = np.delete(v_blocks, [i, j], axis=1)

        # print(b_new)




    return s_blocks[0]