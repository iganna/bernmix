# import numpy as np
# import re
#
#
#
# def int_2_multbern(x):
#     x = int(x)
#     sx = bin(x)[2:]
#     y = re.findall(r'\d', sx)
#     # a = [i for i, y1 in enumerate(y) if y1 == '1']
#
#     indiv =
#
#     return
#
#
#
#
# def conv_population(weights, target, n_sol):
#     """
#     A Dynamic Programming solution for subset sum problem
#     :param weights:
#     :param target:
#     :param n_sol:
#     :return:
#     """
#
#     solutions = np.zeros((target + 1, n_sol))
#     outcomes = np.zeros(target + 1).astype(int)
#     outcomes[0] = 1
#     max_outcome = 1
#
#     for i, w in enumerate(weights):
#         print(i)
#         for j in range(max_outcome):
#             k = w + j
#             if k > target:
#                 break
#             if outcomes[j] == 0:
#                 continue
#             for s in solutions[j, :outcomes[j]]:
#                 if outcomes[k] == n_sol:
#                     break
#                 solutions[k, outcomes[k]] = s + 2 ** i
#                 outcomes[k] += 1
#                 max_outcome = max(max_outcome, k)
#         if outcomes[target] == n_sol:
#             break
#
#
# n = 100
# m = 10**4
# weights = np.random.rand(n)
# weights = weight_rounding(weights, m)
# target_indiv = np.random.binomial(1, 0.1, n)
# print(sum(target_indiv))
# target = np.dot(target_indiv, weights)