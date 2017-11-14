from scipy import optimize as opt
import numpy as np
import timeit

n = 100
w = np.random.rand(n)
M = 10 ** 4
coef = (M + 2*n) / sum(w)
wc = np.around(w * coef).astype(int)

target_indiv = np.random.randint(0, 2, n).astype(int)
target_value_z = np.dot(wc, target_indiv).astype(int)



a = (np.random.rand(n)-0.5)*2



c = np.append(wc, -target_value_z)

A_ident = np.identity(n + 1).astype(int)
A_ident = np.delete(A_ident, n, axis=0)
A_ub = np.concatenate((A_ident, [-c]), axis=0)
b_ub = np.append(np.ones(n), 0).astype(int)


# x_(n+1) = 1
A_eq = [np.append(np.zeros(n), 1).astype(int)]
b_eq = 1
# idx = np.random.permutation(n)
# idx = [15]
# for i in idx:
#     tmp = np.zeros(n+1);
#     tmp[i] = 1
#     A_eq = np.concatenate((A_eq, [tmp]), axis=0)
#     b_eq = np.append(b_eq, np.random.randint(0, 2, 1))
#

x = opt.linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq)
x.x
np.dot(x.x, c)
np.dot(np.around(x.x), c)

# bounds_up = np.ones(n+1).astype(int)
# bounds_down = np.append(np.zeros(n), 1).astype(int)
#
# x = opt.linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=(bounds_down, bounds_up))
#
#
# x = opt.linprog(c, A_ub=[-c], b_ub=0, bounds=(bounds_down, bounds_up))

target_value_s = np.dot(w, target_indiv) * coef
k = 5
n_step = 1000
step = 2*k / n_step
target_range = np.arange(target_value_s-k, target_value_s+k, step)

pop = np.empty((0, n+1))
start = timeit.timeit()
for target_value in target_range:
    c = np.append(w * coef, -target_value)

    A_ident = np.identity(n + 1).astype(int)
    A_ident = np.delete(A_ident, n, axis=0)
    A_ub = np.concatenate((A_ident, [-c]), axis=0)
    b_ub = np.append(np.ones(n), 0).astype(int)

    # x_(n+1) = 1
    A_eq = [np.append(np.zeros(n), 1).astype(int)]
    b_eq = 1
    x = opt.linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq)
    pop = np.concatenate((pop, [x.x]), axis=0)

    np.dot(x.x, c)


end = timeit.timeit()
print(end - start)


