import os
import numpy as np
import daiquiri

import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions
tfb = tfp.bijectors

from scipy.optimize import minimize

from omega import __logger_name__, __version__
logger = daiquiri.getLogger(__logger_name__ + '.estimator.omega')

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 


def get_reparameterized_negative_binomial(mean, overdispersion):
    
    """
    This parameterization is consistent with: 
    variance = mean + overdispersion * mean^2
    reference: https://github.com/tensorflow/probability/issues/372
    """
    
    total_count = 1 / overdispersion
    p = mean / (mean + total_count)
    return tfd.NegativeBinomial(total_count, probs=p)


def get_synthetic_data(l, true_omega, true_overdispersion):

    neg_binom = get_reparameterized_negative_binomial(l * true_omega, true_overdispersion)
    return neg_binom.sample()


def dichotomous_search(point_estimate, func, bound, step=5, tol=1e-3):
    
    # step 1: looking for the upper limit
    a = point_estimate
    b = a + step
    diff = bound - func(b)
    while diff > 0:
        a = b
        b += step
        diff = bound - func(b)
    eps = abs(a - b)
    
    # step 2: refine upper limit with dichotomous search
    while eps > tol:
        middle = (a + b) / 2
        if (bound - func(middle)) * (bound - func(b)) > 0:
            b = middle
        else:
            a = middle
        eps = abs(a - b)
    
    upper_limit = a

    # step 3: set a lower limit
    a = point_estimate
    b = 0.01
    eps = abs(a - b)
    
    # step 4: refine lower limit with dichotomous search
    while eps > tol:
        middle = (a + b) / 2
        if (bound - func(middle)) * (bound - func(b)) > 0:
            b = middle
        else:
            a = middle
        eps = abs(a - b)
    
    lower_limit = a

    return lower_limit, upper_limit

    
@tf.function(autograph=False, experimental_compile=True)
def sampler(num_results, num_burnin_steps, log_prob_func):
    
    hmc = tfp.mcmc.HamiltonianMonteCarlo(
        target_log_prob_fn=log_prob_func,
        num_leapfrog_steps=10,
        step_size=0.05
    )

    initial_state = [
        tf.ones([], name='init_omega')
    ]
    
    unconstraining_bijectors = [
        tfb.Exp()
    ]

    kernel = tfp.mcmc.TransformedTransitionKernel(
        inner_kernel=hmc, 
        bijector=unconstraining_bijectors
    )

    samples = tfp.mcmc.sample_chain(
        num_results=num_results,
        num_burnin_steps=num_burnin_steps,
        current_state=initial_state,
        kernel=kernel,
        trace_fn=None
    )

    return samples


class dNdS:

    def __init__(self, l, n, dispersion):

        # Get the indices of non-zero values in l
        non_zero_indices = tf.where(tf.not_equal(l, 0))

        # Use tf.gather_nd to gather values from n using non-zero indices
        l_non_zero = tf.gather_nd(l, non_zero_indices)
        n_non_zero = tf.gather_nd(n, non_zero_indices)
        
        self.l = l_non_zero
        self.n = n_non_zero

        # only relevant for mle, in bayes defaults to 1 but is not used
        self.dispersion = dispersion
        
        self.res = None
        self.vector_size = len(non_zero_indices)

    def mle_run(self, debug=False):
        
        dispersion = self.dispersion
        
        def minus_log_like(w):
            mu = w * self.l
            if dispersion > 0:
                f = 1 / dispersion
                p = mu / (mu + f)
                model = tfd.NegativeBinomial(f, probs=p)
            elif dispersion == 0:
                model = tfd.Poisson(mu)
            return -tf.reduce_sum(model.log_prob(self.n))

        res = minimize(minus_log_like, 1., method='nelder-mead', options={'xatol': 1e-8, 'disp': False})
        omega_hat = res.x[0]

        def twice_llr(w):
            return 2 * (minus_log_like(w) - minus_log_like(omega_hat))

        # MLE log-likelihood
        l1 = minus_log_like(omega_hat)
        l0 = minus_log_like(1.)
        
        # LRT
        lambda_ = 2 * (l0 - l1)
        pvalue = tfd.Chi2(1.).survival_function(lambda_)

        # Confidence intervals
        alpha = 0.05
        chi2 = tfp.distributions.Chi2(1)
        llr_boundary = chi2.quantile(1-alpha).numpy()
        lower, upper = dichotomous_search(omega_hat, twice_llr, llr_boundary)

        return omega_hat, lower, upper, pvalue.numpy()


    def bayes_run(self, debug=False):
        
        # Lognormal-Poisson model
        model = tfd.JointDistributionSequential([
            # dN/dS prior: https://en.wikipedia.org/wiki/Log-normal_distribution#/media/File:Log-normal-pdfs.png
            # TODO: we want a non-informative prior centered at ~1, how skewed towards >= 1?
            # n: mutation count
            # tfd.LogNormal(loc=0., scale=1.),
            tfd.LogNormal(loc=0., scale=0.25),
            lambda dnds: tfd.Poisson(self.l * dnds),
            ])

        def log_prob_func(omega):
            return tf.reduce_mean(model.log_prob([omega, self.n]))
        
        num_results = 10000
        num_burnin_steps = 100000
        chain = sampler(num_results, num_burnin_steps, log_prob_func)

        return chain


def bayes_infer(args):
    
    gene_term, sample_term, impact_term, gene_set, sample_set, impact_set, l, n = args

    res = {}

    dnds_calculator = dNdS(l, n, 1)


    res['gene'] = [gene_term]
    res['sample'] = [sample_term]
    res['impact'] = [impact_term]
    res['mutations'] = [int(sum(n).numpy())]

    try:
        chain = dnds_calculator.bayes_run()
        res['mean_dnds'] = [np.mean(chain)]
        res['perc_25_dnds'] = [np.percentile(chain, 25)]
        res['perc_75_dnds'] = [np.percentile(chain, 75)]

    except Exception as e:
        logger.warning(f"Bayes method did not work because of {e}")
        res['mean_dnds'] = [None]
        res['perc_25_dnds'] = [None]
        res['perc_75_dnds'] = [None]
    
    return res


def mle_infer(args, dispersion):

    gene_term, sample_term, impact_term, gene_set, sample_set, impact_set, l, n = args

    res = {}
    res_learning_curve = {}

    dnds_calculator = dNdS(l, n, dispersion)
    # logger.debug(f"dNdS calculator for {gene_term}\t{sample_term}\t{impact_term}\n{gene_set}\t{sample_set}\t{impact_set}\nhas n equal to {n} and l equal to {l}\nn, l pairs\n{list(zip(list(n), list(l)))}")

    res['gene'] = [gene_term]
    res['sample'] = [sample_term]
    res['impact'] = [impact_term]
    res['mutations'] = [int(sum(n).numpy())]

    res_learning_curve['gene'] = [gene_term]
    res_learning_curve['sample'] = [sample_term]
    res_learning_curve['impact'] = [impact_term]

    # try:
    omega_hat, lower, upper, pvalue = dnds_calculator.mle_run()

    res['dnds'] = [omega_hat]
    res['pvalue'] = [pvalue]
    res['lower'] = [lower]
    res['upper'] = [upper]

    # TODO: remove res_learning_curve
    res_learning_curve['learning_curve'] = [None]

    # except:
    #     res['dnds'] = [None]
    #     res['pvalue'] = [None]
    #     res['lower'] = [None]
    #     res['upper'] = [None]

    #     res_learning_curve['learning_curve'] = [[None]*1000]

    return res, res_learning_curve


if __name__ == '__main__':
    pass
