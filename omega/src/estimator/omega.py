import os
import numpy as np
import daiquiri

import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions
tfb = tfp.bijectors

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

    def __init__(self, l, n):
        
        self.l = l
        self.n = n
        self.res = None

    def mle_run(self, debug=False):
        
        # dN/dS parameter
        omega = tfp.util.TransformedVariable(1., tfp.bijectors.Exp(), name='omega')

        # instantiate negative binomial model
        dispersion = 0.1
        f = 1 / dispersion
        mean = tfp.util.DeferredTensor(omega, lambda x: self.l * x, shape=(96,))
        p = tfp.util.DeferredTensor(mean, lambda x: x / (x + f), shape=(96,))
        model = tfd.NegativeBinomial(f, probs=p)

        # null log-likelihood
        l0 = -tf.reduce_sum(model.log_prob(self.n))
        
        # regularization parameter
        # alpha = 10  

        # learning rate schedule
        lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
            initial_learning_rate=1e-2,
            decay_steps=1000,
            decay_rate=0.9)

        # convergence criterion
        convergence_criterion = tfp.optimizer.convergence_criteria.LossNotDecreasing(
            rtol=0.1, window_size=1, min_num_steps=25)


        # MLE optimization
        self.res = tfp.math.minimize(
            # loss_fn=lambda: -tf.reduce_sum(model.log_prob(self.n)) + alpha * omega,
            loss_fn=lambda: -tf.reduce_sum(model.log_prob(self.n)),
            num_steps=1000,
            convergence_criterion=convergence_criterion,
            optimizer=tf.optimizers.Adam(learning_rate=lr_schedule),
            trainable_variables=model.trainable_variables
        )

        # MLE omega estimate
        omega_hat = tf.convert_to_tensor(omega)

        def log_like(w):

            f = 1 / dispersion
            mu = w * self.l
            p = mu / (mu + f)
            model = tfd.NegativeBinomial(f, probs=p)
            return tf.reduce_sum(model.log_prob(self.n))

        def twice_llr(w):

            return 2 * (log_like(omega_hat) - log_like(w))

        # MLE log-likelihood
        l1 = -tf.reduce_sum(model.log_prob(self.n))

        # LRT
        lambda_ = 2 * (l0 - l1)
        pvalue = tfd.Chi2(1.).survival_function(lambda_)

        # Confidence intervals

        alpha = 0.05
        chi2 = tfp.distributions.Chi2(1)
        llr_boundary = chi2.quantile(1-alpha).numpy()
        lower, upper = dichotomous_search(omega_hat, twice_llr, llr_boundary)

        return omega_hat.numpy(), lower.numpy(), upper.numpy(), pvalue.numpy(), self.res.numpy()


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
    
    dnds_calculator = dNdS(l, n)

    res['gene'] = [gene_term]
    res['sample'] = [sample_term]
    res['impact'] = [impact_term]

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


def mle_infer(args):

    gene_term, sample_term, impact_term, gene_set, sample_set, impact_set, l, n = args

    res = {}
    res_learning_curve = {}
    
    dnds_calculator = dNdS(l, n)

    res['gene'] = [gene_term]
    res['sample'] = [sample_term]
    res['impact'] = [impact_term]

    res_learning_curve['gene'] = [gene_term]
    res_learning_curve['sample'] = [sample_term]
    res_learning_curve['impact'] = [impact_term]

    # try:
    omega_hat, lower, upper, pvalue, learning_curve = dnds_calculator.mle_run()

    res['dnds'] = [omega_hat]
    res['pvalue'] = [pvalue]
    res['lower'] = [lower]
    res['upper'] = [upper]

    res_learning_curve['learning_curve'] = [learning_curve]

    # except:
    #     res['dnds'] = [None]
    #     res['pvalue'] = [None]
    #     res['lower'] = [None]
    #     res['upper'] = [None]

    #     res_learning_curve['learning_curve'] = [[None]*1000]

    return res, res_learning_curve


if __name__ == '__main__':
    pass
