import os
import itertools

import numpy as np

import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions
tfb = tfp.bijectors

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 


def trinucleotide_contexts():
    
    cb = dict(zip('ACGT', 'TGCA'))
    subs = [''.join(z) for z in itertools.product('CT', 'ACGT') if z[0]!=z[1]]
    flanks = [''.join(z) for z in itertools.product('ACGT', repeat=2)]
    contexts_unformatted = sorted([(a, b) for a, b in itertools.product(subs, flanks)], key=lambda x: (x[0], x[1]))
    contexts = [b[0]+a[0]+b[1]+'>'+a[1] for a, b in contexts_unformatted]
    return contexts


def transform_bracket_context(bracket_context):

    ref = bracket_context[2]
    alt = bracket_context[4]
    flank1 = bracket_context[0]
    flank2 = bracket_context[-1]
    return flank1 + ref + flank2 + '>' + alt


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
        total_count = 1 / dispersion
        mean = tfp.util.DeferredTensor(omega, lambda x: self.l * x, shape=(96,))
        p = tfp.util.DeferredTensor(mean, lambda x: x / (x + total_count), shape=(96,))
        model = tfd.NegativeBinomial(total_count, probs=p)

        # null log-likelihood
        l0 = -tf.reduce_sum(model.log_prob(self.n))
        
        # regularization parameter
        alpha = 10  

        # MLE optimization
        self.res = tfp.math.minimize(
            loss_fn=lambda: -tf.reduce_sum(model.log_prob(self.n)) + alpha * omega,
            num_steps=50,
            optimizer=tf.optimizers.Adam(learning_rate=0.05), 
            trainable_variables=model.trainable_variables
        )
        # MLE omega estimate
        omega_hat = tf.convert_to_tensor(omega)

        # MLE log-likelihood
        l1 = -tf.reduce_sum(model.log_prob(self.n))
        
        # LRT
        lambda_ = 2 * (l0 - l1)
        pvalue = tfd.Chi2(1.).survival_function(lambda_)

        return omega_hat.numpy(), pvalue.numpy()
    
    def bayes_run(self, debug=False):
        
        # Lognormal-Poisson model
        model = tfd.JointDistributionSequential([
            tfd.LogNormal(loc=0., scale=1.),     
            # dN/dS: https://en.wikipedia.org/wiki/Log-normal_distribution#/media/File:Log-normal-pdfs.png
            # TODO: we want a non-informative prior centered at ~1, how skewed towards >= 1?
            # n: mutation count
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
    except:
        res['mean_dnds'] = [None]
        res['perc_25_dnds'] = [None]
        res['perc_75_dnds'] = [None]
    
    return res


def mle_infer(args):

    gene_term, sample_term, impact_term, gene_set, sample_set, impact_set, l, n = args

    res = {}
    
    dnds_calculator = dNdS(l, n)

    res['gene'] = [gene_term]
    res['sample'] = [sample_term]
    res['impact'] = [impact_term]

    try:
        omega_hat, pvalue = dnds_calculator.mle_run()
        res['dnds'] = [omega_hat]
        res['pvalue'] = [pvalue]
    except:
        res['dnds'] = [None]
        res['pvalue'] = [None]
        
    return res


if __name__ == '__main__':
    pass
