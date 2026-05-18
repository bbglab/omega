import os
import daiquiri

from scipy.optimize import minimize

from omega import __logger_name__, __version__

SEED = 31

logger = daiquiri.getLogger(__logger_name__ + '.estimator.omega')
MIN_PVALUE_FLOAT32_RESOLUTION = 1.17e-38

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['TF_DETERMINISTIC_OPS'] = 'true'
os.environ['TF_CUDNN_DETERMINISTIC']= '1'

def get_reparameterized_negative_binomial(mean, overdispersion):
    """
    This parameterization is consistent with:
    variance = mean + overdispersion * mean^2
    reference: https://github.com/tensorflow/probability/issues/372
    """
    import tensorflow_probability.python.distributions as tfd

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



class dNdS:

    def __init__(self, l, n, dispersion):

        self.tf, self.tfd, self.tfb, self.tfp_mcmc = self._lazy_import_tf()

        # Get the indices of non-zero values in l
        non_zero_indices = self.tf.where(self.tf.not_equal(l, 0))

        # Use tf.gather_nd to gather values from n using non-zero indices
        l_non_zero = self.tf.gather_nd(l, non_zero_indices)
        n_non_zero = self.tf.gather_nd(n, non_zero_indices)

        self.l = l_non_zero
        self.n = n_non_zero

        # only relevant for mle, in bayes defaults to 1 but is not used
        self.dispersion = self.tf.Variable(dispersion)

        self.res = None
        self.vector_size = len(non_zero_indices)

    def _lazy_import_tf(self):
        """Lazy import TensorFlow only when needed"""
        import tensorflow as tf
        import tensorflow_probability as tfp

        tfd = tfp.distributions
        tfb = tfp.bijectors
        tfp_mcmc = tfp.mcmc

        tf.random.set_seed(SEED)
        tfp.util.SeedStream(SEED, salt="random_beta")

        return tf, tfd, tfb, tfp_mcmc


    def sampler(self, num_results, num_burnin_steps, log_prob_func):

        hmc = self.tfp_mcmc.HamiltonianMonteCarlo(
            target_log_prob_fn=log_prob_func,
            num_leapfrog_steps=10,
            step_size=0.05
        )

        initial_state = [
            self.tf.ones([], name='init_omega')
        ]

        unconstraining_bijectors = [
            self.tfb.Exp()
        ]

        kernel = self.tfp_mcmc.TransformedTransitionKernel(
            inner_kernel=hmc,
            bijector=unconstraining_bijectors
        )

        samples = self.tfp_mcmc.sample_chain(
            num_results=num_results,
            num_burnin_steps=num_burnin_steps,
            current_state=initial_state,
            kernel=kernel,
            trace_fn=None
        )

        return samples

    def mle_run(self, debug=False):
        dispersion = self.dispersion

        def minus_log_like(w):
            mu = w * self.l
            if dispersion > 0:
                f = 1 / dispersion
                p = mu / (mu + f)
                model = self.tfd.NegativeBinomial(f, probs=p)
            else:
                model = self.tfd.Poisson(mu)
            return -self.tf.reduce_sum(model.log_prob(self.n))

        res = minimize(minus_log_like, 1., method='nelder-mead', options={'xatol': 1e-8, 'disp': False})
        omega_hat = res.x[0]

        def twice_llr(w):
            return 2 * (minus_log_like(w) - minus_log_like(omega_hat))

        # MLE log-likelihood
        l1 = minus_log_like(omega_hat)
        l0 = minus_log_like(1.)

        # LRT
        lambda_ = 2 * (l0 - l1)
        pvalue = self.tfd.Chi2(1.).survival_function(lambda_)

        # Confidence intervals
        alpha = 0.05
        chi2 = self.tfd.Chi2(1)
        llr_boundary = chi2.quantile(1-alpha).numpy()
        lower, upper = dichotomous_search(omega_hat, twice_llr, llr_boundary)

        return omega_hat, lower, upper, pvalue.numpy()


    def bayes_run(self, debug=False):
        
        # Lognormal-Poisson model
        model = self.tfd.JointDistributionSequential([
            # dN/dS prior: https://en.wikipedia.org/wiki/Log-normal_distribution#/media/File:Log-normal-pdfs.png
            # TODO: we want a non-informative prior centered at ~1, how skewed towards >= 1?
            # n: mutation count
            # tfd.LogNormal(loc=0., scale=1.),
            self.tfd.LogNormal(loc=0., scale=0.25),
            lambda dnds: self.tfd.Poisson(self.l * dnds),
            ])

        def log_prob_func(omega):
            return self.tf.reduce_mean(model.log_prob([omega, self.n]))
        
        num_results = 10000
        num_burnin_steps = 100000
        chain = self.sampler(num_results, num_burnin_steps, log_prob_func)

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

    gene_term, sample_term, impact_term, _gene_set, _sample_set, _impact_set, l, n = args

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

    omega_hat, lower, upper, pvalue = dnds_calculator.mle_run()

    res['dnds'] = [omega_hat]
    pvalue_value = float(pvalue)
    if pvalue_value < MIN_PVALUE_FLOAT32_RESOLUTION:
        logger.warning(
            'P-value (%g) below float32 resolution for gene=%s sample=%s impact=%s; replacing with %g',
            pvalue_value,
            gene_term,
            sample_term,
            impact_term,
            MIN_PVALUE_FLOAT32_RESOLUTION,
        )
        pvalue_value = MIN_PVALUE_FLOAT32_RESOLUTION
    res['pvalue'] = [pvalue_value]
    res['lower'] = [lower]
    res['upper'] = [upper]

    return res
