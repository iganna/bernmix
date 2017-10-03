import numpy as np
import itertools as it
import bernmix_double.bernvect as bv

class BernoulliPop:
    """
    This class describes a population of individuals.
    Individuals are a mutable with a mutation rate, mut_rate
    """

    def __init__(self, indiv_len, pop_size, mut_prob=None, pop_init=None):

        def initial_pop(n):
            # ANNA: np.unique
            return np.random.randint(0, 2, n)

        if not isinstance(indiv_len, int) or indiv_len < 0:
            raise ValueError('n, length of individual must be a non-negative integer')
        if not isinstance(pop_size, int) or pop_size < 0:
            raise ValueError('pop_size must be a non-negative integer')

        # the number of Bernoulli random variables representing an individual
        self._indiv_len = indiv_len
        self._pop_size = pop_size

        self._mut_prob = mut_prob or 2/indiv_len
        self._pop = pop_init or initial_pop((self._pop_size, self._indiv_len))

    @property
    def indiv_len(self):
        return self._indiv_len

    @property
    def mut_prob(self):
        return self._mut_prob

    @property
    def pop_size(self):
        return self._pop_size

    @property
    def pop(self):
        return self._pop

    @pop.setter
    def pop(self, new_pop):
        # ANNA: add check
        self._pop = new_pop

    def rand_indiv(self):
        return np.random.randint(0, 2, self.indiv_len)

    def generate(self):
        """
        This function copy and mutate the population
        :return:
        """

        def mutate(new_pop):
            """
            This function mutate the individual
            :param new_pop:
            :return:
            """
            for i, j in it.product(range(self.pop_size), range(self.indiv_len)):
                if np.random.random() <= self.mut_prob:
                    new_pop[i, j] = 1 - new_pop[i, j]
            return new_pop

        add_pop = mutate(np.copy(self.pop))
        self.pop = np.concatenate((self.pop, add_pop), axis=0)
        _, idx = np.unique(self.pop, return_index=True, axis=0)
        self.pop = self.pop[np.sort(idx)]

    def remove(self, indexes):
        self.pop = np.delete(self.pop, indexes, axis=0)


class EvolvingPop(BernoulliPop):
    """
    The population consists of many individuals with the corresponding weighted sum
    """

    def __init__(self, weights, pop_size, mut_prob=None):
        # ANNA: add asserts
        super().__init__(len(weights), pop_size, mut_prob)
        self._weights = weights
        # self._pop_a = None
        # fitness
        self._fitness = None

    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, weights_new):
        self._weights = weights_new

    @property
    def fitness(self):
        return self._fitness

    @fitness.setter
    def fitness(self, new_fitness):
        # ANNA asserts
        self._fitness = new_fitness

    def weighted_sum(self, indiv, weights=None):
        weights = self.weights if weights is None else weights
        return np.dot(weights, indiv)

    def fitness_func(self, target_value, indexes=None):
        indexes = indexes or range(self.pop_size)
        return list(map(lambda x: abs(x-target_value), map(self.weighted_sum, self.pop[indexes])))

    def fitness_func_signed(self, target_value, indexes=None):
        indexes = indexes or range(self.pop_size)
        return list(map(lambda x: x-target_value, map(self.weighted_sum, self.pop[indexes])))

    def evolve(self, n_epoch, target_value):
        self.fitness = np.array(self.fitness_func(target_value))
        # for _ in it.repeat(None, n_epoch):
        for i_epoch in range(n_epoch):

            if i_epoch == round(i_epoch/100)*100:
                print(i_epoch, sum(self.fitness))
            # duplicate and mutate
            self.generate()
            # calculate fitness of new individuals and add
            self.fitness = np.concatenate((self.fitness,
                                           self.fitness_func(target_value, range(self.pop_size, len(self.pop)))))
            # get indexes of individuals with best fitness
            indexes = self.fitness.argsort()[:self.pop_size]
            # remove individuals with low fitness
            self.pop = self.pop[indexes]
            self.fitness = self.fitness[indexes]

    def evolve_new(self, n_epoch, target_value):

        fitness_signed = np.array(self.fitness_func_signed(target_value))
        (self.pop, fitness_signed_new) = bv.evolving(self.pop,
                               fitness_signed,
                               self.weights,
                               self.mut_prob,
                               n_epoch)

        # print('RMSE genetic olgorithm', sum(abs(fitness_signed_new -
        #                                                 np.array(self.fitness_func_signed(target_value)))))
        self.fitness = np.array(self.fitness_func(target_value))


    def split(self, target_value):
        """
        Stplit population into two populations around target_value:
        """
        values = np.array(list(map(self.weighted_sum, self.pop)))
        indexes_down = np.where(values <= target_value)[0]
        indexes_up = np.where(values > target_value)[0]
        return indexes_down, indexes_up

    def approximate(self, m1, m2, target_value):
        """
        This functions findes the optimal parameter m
        to approximate the distribution at the target_value
        :param m1:
        :param m2:
        :param target_value:
        :return:
        """
        def overlap(set_down, set_up):
            return sum(elem >= min(set_up) for elem in set_down) + \
                   sum(elem <= max(set_down) for elem in set_up)
        # # Alternative error function
        # def overlap_dist(set_down, set_up):
        #     return  max(set_down) - min(set_up)


        def error_func(n, w, m):
            # scaling coefficient
            c = (m + n) / w
            # rounding
            weights_new = np.around(self.weights * c, decimals=1).astype(int)

            values_down = list(map(lambda x: self.weighted_sum(x, weights_new),
                                   self.pop[indexes_down]))
            values_up = list(map(lambda x: self.weighted_sum(x, weights_new),
                                 self.pop[indexes_up]))
            n_overlap = overlap(values_down, values_up)
            return m, n_overlap

        indexes_down, indexes_up = self.split(target_value)

        n = self.indiv_len
        w = sum(self.weights)

        # m is the number of point to approximate distribution
        error_params = map(lambda m: error_func(n, w, m), range(m1, m2+1))
        m_optimal, n_overlap = min(error_params, key=lambda x: x[1])

        print('min Number of overlapped points is', m_optimal, n_overlap)

        return m_optimal, indexes_down, indexes_up

    def prob_discrepancy(self,
                         n_points,
                         indexes_down,
                         indexes_up,
                         target_indiv,
                         probs):
        """
        Calculate the probability that should be relaxed
        :param m_optimal:
        :param indexes_down:
        :param indexes_up:
        :param target_value:
        :return:
        """
        def prob_of_indiv(indiv, prob):
            return np.prod(prob[indiv == 1])


        c = (n_points + self.indiv_len) / sum(self.weights)
        # rounding
        weights_new = np.around(self.weights * c, decimals=1).astype(int)
        target_value = self.weighted_sum(target_indiv, weights_new)
        values_down = list(map(lambda x: self.weighted_sum(x, weights_new),
                               self.pop[indexes_down]))
        values_up = list(map(lambda x: self.weighted_sum(x, weights_new),
                             self.pop[indexes_up]))

        idx_down_descr = indexes_down[values_down > target_value]
        idx_up_descr = indexes_up[values_up <= target_value]

        # probabilities of individuals that should be accounted in sum
        prob_down_descr = list(map(lambda indiv: prob_of_indiv(indiv, probs), idx_down_descr))
        # probabilities of individuals that should NOT be accounted in sum
        prob_up_descr = list(map(lambda indiv: prob_of_indiv(indiv, probs), idx_up_descr))

        print(prob_down_descr)
        print(prob_up_descr)
        return sum(prob_down_descr) - sum(prob_up_descr)

# ---------------------------
#         Pipeline
# ---------------------------
if __name__ == "__main__":

    # initialisation
    weights = np.random.rand(20)
    pop_size = 1000
    mut_prob = 0.2
    n_epoch = 100


    pop = EvolvingPop(weights, pop_size, mut_prob)
    target_indiv = pop.rand_indiv()
    target_value = pop.weighted_sum(target_indiv)
    # target_value = np.random.random() * sum(weights)

    pop.evolve(n_epoch, target_value)

    # parameters for approximation
    m2 = 100000
    m1 = m2 - 1000
    c, x, possible_indiv = pop.approximate(m1, m2, target_value)

