
import bernmix_double.bernmix_double as bmd
import bernmix_int.bernmix_int as bmi



def bernmix_pmf_int(probs, weights):
    """
    This function reputrn the vector of probabilities for possible values
    of the weighted sum of Bernoulli random variables
    when weights are integer
    :param probs:
    :param weights:
    :return: The PMF across all possible values
    """
    pass

def bernmix_cdf_double(probs, weights, target_value, n_points = None):
    """
    This function reputrn the vector of probabilities for possible values
    of the weighted sum of Bernoulli random variables
    when weights are double
    :param probs:
    :param weights:
    :param target_value:
    :param n_points: a number of point to approximate distribution
    :return: the CDF in target_value
    """

    # chech whether the target value is a vector of Bernoulli or the double value

    # if n_point is a single number then - use it
    # if n_point is a list of two - use it as a range for search M

    if n_points is None:
        # envoke genetic algorithm
        pass
    else:
        # evoke with rounding to n_points
        pass

def poibinmix_pmf_int(probs, wights):
    """
    Linear combination of Poisson Binomial distributions
    :param probs:
    :param wights:
    :return:
    """
    pass


def poibinmix_cdf_double(probs, wights, target_value, n_points = None):
    """
    Linear combination of Poisson Binomial distributions
    :param probs:
    :param wights:
    :param target_value:
    :param n_points:
    :return:
    """
    pass


def binmix_pmf_int(probs, num_of_trails, wights):
    """
    Linear combination of Poisson Binomial distributions
    :param probs:
    :param num_of_trails:
    :param wights:
    :return:
    """
    pass

def binmix_cdf_double(probs, num_of_trails, wights, target_value, n_points = None):
    """
    Linear combination of Poisson Binomial distributions
    :param probs:
    :param num_of_trails:
    :param wights:
    :param target_value:
    :param n_points:
    :return:
    """
    pass

def radmix_pmf_int(probs, wights):
    """
    Rademacher distribution
    :return:
    """
    pass




if __name__ == "__main__":
    pass