import numpy as np
import cvxopt
from cvxopt import glpk
import bernmix_int.bernmix_int as bm_int



def cdf_corrected(weights, probs, target_indiv, M, n_solutions = 100):
    '''

    :param weights:
    :param probs:
    :param target_value:
    :return:
    '''


    def prob_of_indiv(indiv, prob):
        n = len(indiv)
        prob_multiply = list(map(lambda i: prob[i] if indiv[i] == 1 else (1 - prob[i]),
                                 range(n)))
        return np.prod(prob_multiply)


    def binprog_multisol(weights, target_value, n_solutions, n_fixed = 10):
        '''
        This function returns multiple solutions of zero-one linear programming problem
        :param c_init:
        :param target_value:
        :param n_solutions:
        :return:
        '''

        glpk.options['msg_lev'] = 'GLP_MSG_OFF'

        n = len(weights)
        c_init = np.append(weights, -target_value).astype(int)
        c = c_init.reshape((1, n+1))

        A_ident = np.identity(n + 1).astype(int)
        A_ident = np.delete(A_ident, n, axis=0)
        A_ub = np.concatenate((A_ident, -c), axis=0)
        b_ub = np.append(np.ones(n), 0).astype(int).reshape(n+1, 1)

        pop = np.empty((0, n+1))
        for _ in range(n_solutions):
            # x_(n+1) = 1
            A_eq = np.append(np.zeros(n), 1).astype(int).reshape(1,n+1)
            b_eq = np.asmatrix(1)

            idx = np.random.permutation(n)
            idx = idx[0:round(n_fixed/100 * n)]
            for i in idx:
                tmp = np.zeros(n+1)
                tmp[i] = 1
                A_eq = np.concatenate((A_eq, [tmp]), axis=0)
                b_eq = np.concatenate((b_eq, [np.random.randint(0, 2, 1)]), axis=0)


            C = cvxopt.matrix(c, tc='d')
            G = cvxopt.matrix(A_ub, tc='d')
            h = cvxopt.matrix(b_ub, tc='d')
            A = cvxopt.matrix(A_eq, tc='d')
            b = cvxopt.matrix(b_eq, tc='d')


            status, x = glpk.ilp(c=C, G=G, h=h, A=A, b=b, B=set(range(n+1)))
            if status != 'optimal':
                continue

            a = np.array(list(xi for xi in x))

            pop = np.concatenate((pop, np.matrix(x).reshape((1, n+1))), axis=0)

        if pop.shape[1] == 0:
            return None

        pop = np.unique(pop, axis=0)
        # ANNA: from 0 to n
        return pop


    n = len(weights)
    w = weights
    coef = (M + 2 * n) / sum(w)
    wc = np.around(w * coef).astype(int)
    target_value_s = np.dot(w, target_indiv)
    target_value_z = np.dot(wc, target_indiv).astype(int)

    # number of points into two directions
    k = 10
    target_range = np.arange(target_value_z - k, target_value_z + k + 1, 1)
    distr_in_probs = np.zeros((2 * k + 1, 2))
    for i, target_value in enumerate(target_range):
        pop = binprog_multisol(wc, target_value, n_solutions, 10)
        if pop is None:
            distr_in_probs[i, 0] = 0
            distr_in_probs[i, 1] = 0
            continue

        s_values = list(map(lambda x: np.dot(x[0:n], w), pop))

        indiv_probs = list(map(lambda x: prob_of_indiv(x[0:n], probs), pop))
        sum_prob = sum(indiv_probs)

        distr_in_probs[i, 0] = sum(np.extract(s_values <= target_value_s, indiv_probs)) / sum_prob
        distr_in_probs[i, 1] = sum(np.extract(s_values > target_value_s, indiv_probs)) / sum_prob

    pmf = bm_int.pmf(probs, wc)

    pmf_range = list(map(lambda t: pmf[t], target_range))
    cdf2 = sum(pmf[:target_value_z + 1 - k - 1:]) + np.dot(pmf_range, distr_in_probs[:, 0])

    return cdf2


if __name__ == "__main__":
    print('bernmix_correction is loaded')