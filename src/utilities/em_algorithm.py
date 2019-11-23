from scipy import stats
import numpy as np
from math import log10, log, isclose
from decimal import Decimal
import matplotlib.pyplot as pyplt

def log_expon_pdf(data, rate):
    """ Pdf for an exponential distribution """
    coeff = float(rate)

    if coeff <= 0:
        raise ValueError("The rate parameter has to be positive.")

    prob = -(data * coeff) + log(coeff)

    return prob


def expon_pdf(data, rate):
    """ Pdf for an exponential distribution """
    coeff = float(rate)

    if coeff <= 0:
        raise ValueError("The rate parameter has to be positive.")

    prob = np.exp(data * (-coeff)) * coeff

    return prob


def estimate_zetas(exponential_likelihood_weighted, uniform_likelihood_weighted):
     zetas = \
         exponential_likelihood_weighted / np.add(exponential_likelihood_weighted,
                                                                     uniform_likelihood_weighted)
     return zetas


def compute_marginal_log_likelihood(exponential_likelihood_weighted, uniform_likelihood_weighted):
    likelihood_sum = np.add(exponential_likelihood_weighted, uniform_likelihood_weighted)
    min_nonzero = np.min(likelihood_sum[np.nonzero(likelihood_sum)])
    likelihood_sum[likelihood_sum == 0] = min_nonzero
    return np.sum(np.log(np.add(exponential_likelihood_weighted, uniform_likelihood_weighted)))


def maximise(x, zetas):
    sum_zeta = float(np.sum(zetas))
    new_phi = sum_zeta / len(zetas)
    new_lambd = sum_zeta / np.sum(zetas * x)

    return new_lambd, new_phi


def em(x, lambd, phi, half_genome_size, max_iter=1000, epsilon=1e-4):
    """ EM algorithm to fix a mixture of an exponential with rate parameter lambd
    and a uniform distribution [0, half_genome_size], with mixture weight phi. """
    uniform_likelihood_single_val = stats.uniform.pdf(x[0], loc=0, scale=half_genome_size)

    marginal_log_likelihoods = []
    i = 0
    has_MLL_not_improved = False

    while (i < max_iter) and has_MLL_not_improved is False:
        exponential_likelihood_weighted = np.exp(log_expon_pdf(x, lambd) + log(phi))
        np.nan_to_num(exponential_likelihood_weighted, copy=False)
        minval = np.min(exponential_likelihood_weighted[np.nonzero(exponential_likelihood_weighted)])
        exponential_likelihood_weighted[exponential_likelihood_weighted == 0] = minval

        uniform_likelihood_weighted = np.repeat(uniform_likelihood_single_val, len(x)) * (1 - phi)

        # Saves the values for this iteration
        marginal_log_likelihoods.append(compute_marginal_log_likelihood(exponential_likelihood_weighted,
                                                                        uniform_likelihood_weighted))

        # Convergency checks
        if i > 0:
            # Prevent numerical instability issues
            has_MLL_not_improved = isclose(marginal_log_likelihoods[i],
                                               marginal_log_likelihoods[i - 1],
                                               rel_tol=epsilon)

        # Estimates zeta
        zetas = estimate_zetas(exponential_likelihood_weighted, uniform_likelihood_weighted)

        # Maximises and updates parameters
        if has_MLL_not_improved == False:
            lambd, phi = maximise(x, zetas)

        i = i + 1


    return {'lambda': lambd,
            'phi': phi,
            'MLL': marginal_log_likelihoods,
            'iteration': i}


def log_sample_parameters(start_value, final_value, n_values):
    if start_value <= 0 or final_value <= 0:
      raise ValueError("The initial and final values have to be greater than zero.")
    if final_value <= start_value:
      raise ValueError("The final value has to be greater than the initial one.")

    log_min = log10(start_value)
    log_max = log10(final_value)
    c = Decimal(log_min)
    i = 0
    while i < n_values:
      yield float(10**c)
      c = c + Decimal((log_max - log_min) / (n_values - 1))
      i = i + 1




def plot_likelihood(likelihood, plot_file_name=None):
    """Function to plot the marginal log likelihood values vs the number
    of iterations and (optional) save it in the file plotFileName."""
    fig = pyplt.figure()
    ax = fig.add_subplot(111)

    ax.plot(likelihood, "o")
    ax.set_yticklabels(ax.get_yticks())
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Marginal Log Likelihood")
    pyplt.tight_layout()

    if plot_file_name is not None:
        pyplt.savefig(plot_file_name)
        pyplt.close()
    else:
        pyplt.show()


def compute_zetas(distances,
                  half_genome_size,
                  lambd,
                  phi):
    uniform_likelihood_single_val = stats.uniform.pdf(distances[0], loc=0, scale=half_genome_size)
    uniform_likelihood_weighted = np.repeat(uniform_likelihood_single_val, len(distances)) * (1 - phi)
    exponential_likelihood_weighted = np.exp(log_expon_pdf(distances, lambd) + log(phi))
    minval = np.min(exponential_likelihood_weighted[np.nonzero(exponential_likelihood_weighted)])
    exponential_likelihood_weighted[exponential_likelihood_weighted == 0] = minval

    return estimate_zetas(exponential_likelihood_weighted, uniform_likelihood_weighted)