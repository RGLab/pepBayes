/*
 * PairedMcmc.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: gimholte
 */

#include "RcppArmadillo.h"
#include "RngStream.h"
#include "Rv.h"
#include <omp.h>
#include "Densities.h"
#include "AdaptiveGibbsSampler.h"
#include "FitBase.h"
#include "PairedMcmc.h"

PairedMcmc::PairedMcmc(const Rcpp::List & par_list, const Rcpp::List & data_list,
        const Rcpp::List & chain_pars) : MarkovChain(chain_pars),
        nu_beta0_shape(Rcpp::as<double>(par_list["nu_beta0_shape"])),
        nu_beta0_rate(Rcpp::as<double>(par_list["nu_beta0_rate"])),
        m_beta0_prior_mean(Rcpp::as<double>(par_list["m_beta0_prior_mean"])),
        m_beta0_prior_prec(Rcpp::as<double>(par_list["m_beta0_prior_prec"])),
        nu_beta1_shape(Rcpp::as<double>(par_list["nu_beta1_shape"])),
        nu_beta1_rate(Rcpp::as<double>(par_list["nu_beta1_rate"])),
        m_beta1_prior_mean(Rcpp::as<double>(par_list["m_beta1_prior_mean"])),
        m_beta1_prior_prec(Rcpp::as<double>(par_list["m_beta1_prior_prec"])),
        lambda_alpha_rate(Rcpp::as<double>(par_list["lambda_alpha_rate"])),
        lambda_eps_rate(Rcpp::as<double>(par_list["lambda_eps_rate"])),
        s_alpha_rate(Rcpp::as<double>(par_list["s_alpha_rate"])),
        s_eps_rate(Rcpp::as<double>(par_list["s_eps_rate"])),
        lambda_a(Rcpp::as<double>(par_list["lambda_a"])),
        lambda_b(Rcpp::as<double>(par_list["lambda_b"])),
        m_beta1_tuner(n_burn),
        nu_beta1_tuner(n_burn)
{
    // set seed for arma's RNG generator so that initializations
    // are reproducible.
    std::srand(seed);

    y_mean = Rcpp::as<arma::mat>(data_list["y_mean"]);
    y_ss = Rcpp::as<arma::mat>(data_list["y_ss"]);
    n_rep = Rcpp::as<arma::ivec>(data_list["n_rep"]);
    pos = Rcpp::as<arma::ivec>(data_list["pos"]);
    pos_start = Rcpp::as<arma::ivec>(data_list["pos_start"]);
    pos_count = Rcpp::as<arma::ivec>(data_list["pos_count"]);

    n_position = pos_start.n_elem;
    n_peptide = y_ss.n_rows;
    n_slide = y_ss.n_cols;
    n_subject = n_slide / 2;

    beta0 = Rcpp::as<arma::vec>(par_list["beta0"]);
    beta1 = Rcpp::as<arma::vec>(par_list["beta1"]);
    mu = Rcpp::as<arma::vec>(par_list["mu"]);
    alpha0 = Rcpp::as<arma::mat>(par_list["alpha0"]);
    alpha1 = Rcpp::as<arma::mat>(par_list["alpha1"]);

    gamma = Rcpp::as<arma::mat>(par_list["gamma"]);
    gamma_prob = arma::mat(n_peptide, n_subject);
    omega = Rcpp::as<arma::vec>(par_list["omega"]);
    a = Rcpp::as<arma::vec>(par_list["a"]);
    b = Rcpp::as<arma::vec>(par_list["b"]);

    nu_alpha0 = Rcpp::as<arma::vec>(par_list["nu_alpha0"]);
    nu_alpha1 = Rcpp::as<arma::vec>(par_list["nu_alpha1"]);
    nu_eps = Rcpp::as<arma::vec>(par_list["nu_eps"]);

    nu_beta0 = Rcpp::as<double>(par_list["nu_beta0"]);
    m_beta0 = Rcpp::as<double>(par_list["m_beta0"]);
    nu_beta1 = Rcpp::as<double>(par_list["nu_beta1"]);
    m_beta1 = Rcpp::as<double>(par_list["m_beta1"]);
    s_alpha = Rcpp::as<double>(par_list["s_alpha"]);
    lambda_alpha = Rcpp::as<double>(par_list["lambda_alpha"]);
    s_eps = Rcpp::as<double>(par_list["s_eps"]);
    lambda_eps = Rcpp::as<double>(par_list["lambda_eps"]);

    nu_err = Rcpp::as<double>(par_list["nu"]);
    nu_re = nu_err;
    w = Rcpp::as<arma::mat>(par_list["w"]);
    u = Rcpp::as<arma::mat>(par_list["u"]);
    mahala_dist = arma::mat(n_peptide, n_slide, arma::fill::zeros);

    // misc helpers for sampling
    mu_mean = arma::vec(n_slide, arma::fill::zeros);
    mu_temp = arma::vec(n_slide - 1, arma::fill::zeros);
    mu_star_mean = arma::vec(n_slide - 1, arma::fill::zeros);
    mu_star = arma::vec(n_slide - 1, arma::fill::zeros);

    Q = Rcpp::as<arma::mat>(par_list["Q"]);
    Qt = Q.t();
    mu_omega = arma::mat(n_slide - 1, n_slide - 1, arma::fill::zeros);
    mu_prec_inner_diag = arma::vec(n_slide, arma::fill::zeros);

    a_sampler = std::vector<AdaptiveGibbsSampler<BetaHyperConditional> >(n_position);
    b_sampler = std::vector<AdaptiveGibbsSampler<BetaHyperConditional> >(n_position);
    s_alpha_sampler = AdaptiveGibbsSampler<ShapeConditional>(1.0, 20.0);
    s_eps_sampler = AdaptiveGibbsSampler<ShapeConditional>(1.0, 5.0);
    for(int q = 0; q < n_position; q++) {
        a_sampler[q] = AdaptiveGibbsSampler<BetaHyperConditional>(.1, 1.0);
        b_sampler[q] = AdaptiveGibbsSampler<BetaHyperConditional>(1.0, 10.0);
    }

    rng = std::vector<RngStream>(n_threads);
    for(int k = 0; k < n_threads; k++) {
        rng[k] = RngStream_CreateStream("");
    }

    // trace storage
    beta0_trace = arma::mat(n_peptide, n_samples);
    beta1_trace = arma::mat(n_peptide, n_samples);
    mu_trace = arma::mat(n_slide, n_samples);

    omega_trace = arma::mat(n_peptide, n_samples);
    a_trace = arma::mat(n_position, n_samples);
    b_trace = arma::mat(n_position, n_samples);

    nu_eps_trace = arma::mat(n_peptide, n_samples);
    nu_alpha0_trace = arma::mat(n_peptide, n_samples);
    nu_alpha1_trace = arma::mat(n_peptide, n_samples);

    hypers_trace = arma::mat(8, n_samples);

    update_mu = Rcpp::as<bool>(chain_pars["update_mu"]);
    update_beta0 = Rcpp::as<bool>(chain_pars["update_beta0"]);
    update_beta1 = Rcpp::as<bool>(chain_pars["update_beta1"]);
    update_beta_hypers = Rcpp::as<bool>(chain_pars["update_beta_hypers"]);
    update_precision = Rcpp::as<bool>(chain_pars["update_precision"]);
    update_precision_hypers = Rcpp::as<bool>(chain_pars["update_precision_hypers"]);
    update_omega = Rcpp::as<bool>(chain_pars["update_omega"]);
    update_omega_hypers = Rcpp::as<bool>(chain_pars["update_omega_hypers"]);
    update_alpha = Rcpp::as<bool>(chain_pars["update_alpha"]);
    update_weights = Rcpp::as<bool>(chain_pars["update_weights"]);
    update_gamma = Rcpp::as<bool>(chain_pars["update_gamma"]);
    share_nu_alpha = false;//Rcpp::as<bool>(chain_pars["share_nu_alpha"]);

    ppa = arma::mat(n_peptide, n_subject, arma::fill::zeros);
    fitted_response = arma::mat(n_peptide, n_subject, arma::fill::zeros);

    Rcpp::Rcout << "Initialization complete. Running " <<
            (n_samples * n_thin + n_burn) << " total iterations" << std::endl;
    Rcpp::Rcout << "for " << n_subject << " subjects with " << n_peptide <<
            " peptides." << std::endl;
}


void PairedMcmc::iterate() {
    int th_id, p, i, q;
    // the order in which updates occur is important
    #pragma omp parallel private(th_id, p, i, q) num_threads(n_threads)
    {
        th_id = omp_get_thread_num();

        #pragma omp for
        for (p = 0; p < n_peptide; p++) {
            updateGammaAlpha(p, rng[th_id]);
        }

        #pragma omp for
        for (i = 0; i < n_subject; i++) {
            computeMahalaDist(i);
            if (update_weights)
                updateWeights(i, rng[th_id]);
        }

        #pragma omp single
        {
            if (update_beta_hypers)
                updateM0(rng[th_id]);
        } // implied barrier

        #pragma omp for
        for (p = 0; p < n_peptide; p++) {
            if (update_beta0)
                updateBeta0(p, rng[th_id]);
        }

        #pragma omp for
        for (p = 0; p < n_peptide; p++) {
            if (update_beta1)
                updateBeta1(p, rng[th_id]);
        }

        #pragma omp for
        for (i = 0; i < n_subject; i++) {
            if (update_mu)
                computeMuMean(i);
        }

        #pragma omp sections
        {

            #pragma omp section
            {
                if (update_precision_hypers)
                    updateAlphaHypers(rng[th_id]);
            }
            #pragma omp section
            {
                if (update_precision_hypers)
                    updateEpsilonHypers(rng[th_id]);
            }
            #pragma omp section
            {
                if (update_beta_hypers)
                    updateNuBeta0(rng[th_id]);
            }
            #pragma omp section
            {
                if (update_mu)
                    updateMu(rng[th_id]);
            }
            #pragma omp section
            {
                if (update_beta_hypers) {
                    updateM1(rng[th_id]);
                    updateNuBeta1(rng[th_id]);
                }
            }
        }

        #pragma omp for
        for (p = 0; p < n_peptide; p++) {
            if (update_precision)
                updateNuEps(p, rng[th_id]);
        }

        #pragma omp for
        for (p = 0; p < n_peptide; p++) {
            if (update_precision)
                updateNuAlpha(p, rng[th_id]);
        }

        #pragma omp for
        for (q = 0; q < n_position; q++) {
            if (update_omega_hypers)
                updateBetaHypers(q, rng[th_id]);
        }

        #pragma omp for
        for (p = 0; p < n_peptide; p++) {
            if (update_omega)
                updateOmega(p, rng[th_id]);
        }
    }
}

void PairedMcmc::collectIteration(const int & sample_idx) {
    beta0_trace.col(sample_idx) = beta0;
    beta1_trace.col(sample_idx) = beta1;
    mu_trace.col(sample_idx) = mu;
    nu_eps_trace.col(sample_idx) = nu_eps;
    nu_alpha0_trace.col(sample_idx) = nu_alpha0;
    nu_alpha1_trace.col(sample_idx) = nu_alpha1;
    omega_trace.col(sample_idx) = omega;
    a_trace.col(sample_idx) = a;
    b_trace.col(sample_idx) = b;

    hypers_trace(0, sample_idx) = s_alpha;
    hypers_trace(1, sample_idx) = lambda_alpha;
    hypers_trace(2, sample_idx) = s_eps;
    hypers_trace(3, sample_idx) = lambda_eps;
    hypers_trace(4, sample_idx) = m_beta0;
    hypers_trace(5, sample_idx) = nu_beta0;
    hypers_trace(6, sample_idx) = m_beta1;
    hypers_trace(7, sample_idx) = nu_beta1;

    ppa = ppa * ((double) sample_idx / (sample_idx + 1.0)) +
            gamma / (sample_idx + 1.0);
    fitted_response = fitted_response * ((double) sample_idx / (sample_idx + 1.0)) +
            (alpha0 % (1 - gamma) + alpha1) / (sample_idx + 1.0);
    for(unsigned i = 0; i < fitted_response.n_cols; i++) {
        fitted_response.col(i) += (beta0 + beta1 % gamma.col(i)) / (sample_idx + 1.0);
    }
}

const Rcpp::List PairedMcmc::chainOutput() {
    Rcpp::List out = Rcpp::List::create(
            Rcpp::Named("hypers") = hypers_trace,
            Rcpp::Named("beta0") = beta0_trace,
            Rcpp::Named("beta1") = beta1_trace,
            Rcpp::Named("mu") = mu_trace,
            Rcpp::Named("nu_alpha0") = nu_alpha0_trace,
            Rcpp::Named("nu_alpha1") = nu_alpha1_trace,
            Rcpp::Named("nu_eps") = nu_eps_trace,
            Rcpp::Named("omega") = omega_trace,
            Rcpp::Named("a") = a_trace,
            Rcpp::Named("b") = b_trace,
            Rcpp::Named("ppb") = ppa,
            // the remaining parameters are provided for chain restarts
            Rcpp::Named("w") = w,
            Rcpp::Named("u") = u,
            Rcpp::Named("alpha0") = alpha0,
            Rcpp::Named("alpha1") = alpha1,
            Rcpp::Named("gamma") = gamma);
    return out;
}

void PairedMcmc::computeMuMean(const int & i) {
    mu_mean(2 * i) = arma::sum(
            (y_mean.col(2 * i) - beta0 - alpha0.col(i)) % nu_eps % n_rep % w.col(2 * i));
    mu_mean(2 * i + 1) = arma::sum((y_mean.col(2 * i + 1) - beta0 -
            alpha0.col(i) + alpha0.col(i) % gamma.col(i) -
            gamma.col(i) % (alpha1.col(i) + beta1)) % nu_eps % n_rep % w.col(2 * i + 1));
    mu_prec_inner_diag(2 * i) = arma::sum(nu_eps % n_rep % w.col(2 * i));
    mu_prec_inner_diag(2 * i + 1) = arma::sum(nu_eps % n_rep % w.col(2 * i + 1));
}

void PairedMcmc::updateMu(RngStream & rng) {
    // Get cholesky decomposition of conditional precision
    mu_omega = Qt * arma::diagmat(mu_prec_inner_diag) * Q;
    mu_omega.diag() += 1.0;
    mu_omega = arma::chol(mu_omega);
    // do forward/back solve to get the conditional mean
    mu_temp = arma::solve(mu_omega.t(), Qt * mu_mean);
    mu_star_mean = arma::solve(mu_omega, mu_temp);
    // fill with n(0, 1) values
    for(int i = 0; i < n_slide - 1; i++) {
        mu_temp(i) = RngStream_N01(rng);
    }
    // solve linear system to get N(0, omega^(-1)) multivariate.
    mu_star = arma::solve(mu_omega, mu_temp) + mu_star_mean;
    // transform mu_star to recover mu
    mu = Q * mu_star;
    return;
}

void PairedMcmc::updateM0(RngStream & rng) {
    double c_prec(0.0), c_mean(0.0);
    c_mean = arma::sum(beta0) * nu_beta0 +
            m_beta0_prior_mean * m_beta0_prior_prec;
    c_prec = n_peptide * nu_beta0 + m_beta0_prior_prec;
    c_mean /= c_prec;
    m_beta0 = c_mean + RngStream_N01(rng) / sqrt(c_prec);
    return;
}

void PairedMcmc::updateBeta0(const int & p, RngStream & rng) {
    double c_mean(0.0), c_prec(0.0);
    int i;
    for(i = 0; i < n_subject; i++) {
        c_mean += (y_mean(p, 2 * i) - mu(2 * i) - alpha0(p, i)) * w(p, 2 * i) *
                nu_eps(p) * n_rep(p);
        c_mean += (y_mean(p, 2 * i + 1) - mu(2 * i + 1) - alpha0(p, i) -
                gamma(p, i) * (beta1(p) + alpha1(p, i) - alpha0(p, i))) *
                        w(p, 2 * i + 1) * nu_eps(p) * n_rep(p);
        c_prec += w(p, 2 * i) * n_rep(p) * nu_eps(p) +
                w(p, 2 * i + 1) * n_rep(p) * nu_eps(p);
    }
    c_mean += m_beta0 * nu_beta0;
    c_prec += nu_beta0;
    c_mean /= c_prec;

    beta0(p) = c_mean + RngStream_N01(rng) / sqrt(c_prec);
    return;
}

void PairedMcmc::updateBeta1(const int & p, RngStream & rng) {
    double s(0.0), m(0.0), tmp;
    for(int i = 0; i < n_subject; i++) {
        if (gamma(p, i) == 1.0) {
            tmp = w(p, 2 * i + 1);
            m += tmp * (y_mean(p, 2 * i + 1) - mu(2 * i + 1) - beta0(p) -
                    alpha1(p, i));
            s += tmp;
        }
    }
    m *= nu_eps(p) * n_rep(p);
    s *= nu_eps(p) * n_rep(p);
    m += m_beta1 * nu_beta1;
    s += nu_beta1;

    m /= s;
    beta1(p) = RngStream_TruncNorm(m, 1.0 / s, rng);
}

void PairedMcmc::updateAlphaHypers(RngStream & rng) {
    int p;
    double density_args[4];
    double log_sum(0.0);
    double nu_sum(0.0);
    for(p = 0; p < n_peptide; p++) {
        nu_sum += nu_alpha0(p);
        log_sum += log(nu_alpha0(p));
        if (!share_nu_alpha) {
            nu_sum += nu_alpha1(p);
            log_sum += log(nu_alpha1(p));
        }
    }
    double r = share_nu_alpha ? 1.0 * n_peptide : (2.0 * n_peptide);
    density_args[0] = r;
    density_args[1] = log(lambda_alpha_rate + nu_sum);
    density_args[2] = log_sum;
    density_args[3] = s_alpha_rate;
    if(!R_FINITE(log_sum))
        Rcpp::stop("invalid input to alpha hyper sampler");
    s_alpha = s_alpha_sampler.sample(rng, density_args);
    lambda_alpha = RngStream_GA1(1.0 + r * s_alpha, rng) /
            (lambda_alpha_rate + nu_sum);
}

void PairedMcmc::updateEpsilonHypers(RngStream & rng) {
    int p;
    double density_args[4];
    double log_sum(0.0);
    double nu_sum(0.0);
    for(p = 0; p < n_peptide; p++) {
        nu_sum += nu_eps(p);
        log_sum += log(nu_eps(p));
    }
    density_args[0] = (double) n_peptide;
    density_args[1] = log(lambda_eps_rate + nu_sum);
    density_args[2] = log_sum;
    density_args[3] = s_eps_rate;
    if(!R_FINITE(log_sum))
        Rcpp::stop("invalid input to epsilon hyper sampler");
    s_eps = s_eps_sampler.sample(rng, density_args);
    lambda_eps = RngStream_GA1(1.0 + n_peptide * s_eps, rng) /
            (lambda_eps_rate + nu_sum);
}

void PairedMcmc::updateNuBeta0(RngStream & rng) {
    int p;
    double s(0.0);
    for(p = 0; p < n_peptide; p++) {
        s += (beta0(p) - m_beta0) * (beta0(p) - m_beta0);
    }
    s /= 2.0;
    s += nu_beta0_rate;
    nu_beta0 = RngStream_GA1(nu_beta0_shape + n_peptide / 2.0, rng) / s;
}


void PairedMcmc::updateNuEps(const int & p, RngStream & rng) {
    double d0, d1, s(0.0);
    const double np(n_rep(p));
    for(int i = 0; i < n_subject; i++) {
        d0 = y_mean(p, 2 * i) - beta0(p) - mu(2 * i) - alpha0(p, i);
        d1 = y_mean(p, 2 * i + 1) - beta0(p) - mu(2 * i + 1) -
                alpha0(p, i) - gamma(p, i) * (beta1(p) + alpha1(p, i) - alpha0(p, i));
        d0 *= d0 * np;
        d1 *= d1 * np;
        d0 += y_ss(p, 2 * i);
        d1 += y_ss(p, 2 * i + 1);
        d0 *= w(p, 2 * i);
        d1 *= w(p, 2 * i + 1);
        s += d0 + d1;
    }
    s /= 2.0;
    s += lambda_eps;
    nu_eps(p) = RngStream_GA1(s_eps + n_subject * np, rng) / s;
}

void PairedMcmc::updateNuAlpha(const int & p, RngStream & rng) {
    double s0(0.0), s1(0.0), a0, a1, ell(0.0);
    int i;
    for(i = 0; i < n_subject; i++) {
        a0 = alpha0(p, i);
        a1 = alpha1(p, i);
        s0 += u(p, 2 * i) * a0 * a0;
        s1 += u(p, 2 * i + 1) * a1 * a1 * gamma(p, i);
        ell += gamma(p, i);
    }

    if (share_nu_alpha) {
        double l = (s0 + s1) / 2.0 + lambda_alpha;
        double s = s_alpha + (n_subject + ell) / 2.0;
        double nu_tmp = RngStream_GA1(s, rng) / l;
        nu_alpha0(p) = nu_tmp;
        nu_alpha1(p) = nu_tmp;
    } else {
        s0 /= 2.0;
        s0 += lambda_alpha;
        s1 /= 2.0;
        s1 += lambda_alpha;
        nu_alpha0(p) = RngStream_GA1(s_alpha + n_subject / 2.0, rng) / s0;
        nu_alpha1(p) = RngStream_GA1(s_alpha + ell / 2.0, rng) / s1;
    }
}

void PairedMcmc::updateOmega(const int & p, RngStream & rng) {
    const double s = arma::sum(gamma.row(p));
    omega(p) = RngStream_LogitBeta(a(pos(p)) + s, b(pos(p)) + n_subject - s, rng);
}

void PairedMcmc::updateM1(RngStream & rng) {
    const double sum_b1 = arma::sum(beta1);
    const double cond_tau = n_peptide * nu_beta1 + m_beta1_prior_prec;
    const double cond_mu = nu_beta1 * sum_b1 +
            m_beta1_prior_mean * m_beta1_prior_prec;
    const double prop_m1 = m_beta1 * exp((RngStream_RandU01(rng) - .5) *
            m_beta1_tuner.getScale());

    const double cur_lik = - m_beta1 * m_beta1 * cond_tau / 2 +
            cond_mu * m_beta1 -
            n_peptide * R::pnorm(m_beta1 * sqrt(nu_beta1), 0.0, 1.0, 1, 1);
    const double prop_lik = - prop_m1 * prop_m1 * cond_tau / 2 +
            cond_mu * prop_m1 -
            n_peptide * R::pnorm(prop_m1 * sqrt(nu_beta1), 0.0, 1.0, 1, 1);

    const double cur_to_prop = - log(prop_m1);
    const double prop_to_cur = - log(m_beta1);

    m_beta1 = metropolisHastings(m_beta1, cur_lik, cur_to_prop,
            prop_m1, prop_lik, prop_to_cur, rng);
    m_beta1_tuner.update(m_beta1 == prop_m1);
    return;
}

void PairedMcmc::updateNuBeta1(RngStream & rng) {
    const double ss_beta1 = pow(arma::norm(beta1 - m_beta1, 2), 2);
    const double cond_shape = n_peptide / 2.0 + nu_beta1_shape;
    const double cond_rate = ss_beta1 / 2.0 + nu_beta1_rate;
    const double prop_nb1 = nu_beta1 * exp((RngStream_RandU01(rng) - .5) *
            nu_beta1_tuner.getScale());

    const double cur_lik = - nu_beta1 * cond_rate + cond_shape * log(nu_beta1) -
            n_peptide * R::pnorm(m_beta1 * sqrt(nu_beta1), 0.0, 1.0, 1, 1);
    const double prop_lik = - prop_nb1 * cond_rate + cond_shape * log(prop_nb1) -
            n_peptide * R::pnorm(m_beta1 * sqrt(prop_nb1), 0.0, 1.0, 1, 1);

    const double cur_to_prop = -log(prop_nb1);
    const double prop_to_cur = -log(nu_beta1);

    nu_beta1 = metropolisHastings(nu_beta1, cur_lik, cur_to_prop, prop_nb1,
            prop_lik, prop_to_cur, rng);
    nu_beta1_tuner.update(nu_beta1 == prop_nb1);
}

void PairedMcmc::updateBetaHypers(const int & q, RngStream & rng) {
    double a_log_sum = lambda_a;
    double b_log_sum = lambda_b;
    const int p_num = pos_count(q);
    for(int p = 0; p < p_num; p++) {
        a_log_sum -= logFromLogit(omega(p + pos_start(q)));
        b_log_sum -= log1mFromLogit(omega(p + pos_start(q)));
    }
    double density_args[3];
    density_args[0] = (double) p_num;
    density_args[1] = b(q);
    density_args[2] = a_log_sum;
    a(q) = a_sampler[q].sample(rng, density_args);

    density_args[1] = a(q);
    density_args[2] = b_log_sum;
    b(q) = b_sampler[q].sample(rng, density_args);

    if (!R_FINITE(a(q)) || !R_FINITE(b(q)) || a(q) <= 0.0 || b(q) <= 0.0) {
        Rcpp::Rcout << a(q) << " " << b(q) << std::endl;
    }
}

void PairedMcmc::updateGammaAlpha(const int & p, RngStream & rng) {
    // metropolis-hastings update for alpha_{ip1} with an
    // "always-move" proposal
    int i, j;
    double cur_den, prop_den;
    double cur_to_prop, prop_to_cur;
    double cur_val, prop_val, tmp;
    double log_omega = logFromLogit(omega(p));
    double log1m_omega = log1mFromLogit(omega(p));
    double nu_star, alpha_star;
    double lambda_star;

    for (i = 0; i < n_subject; i++) {
        j = 2 * i + 1;
        tmp = y_mean(p, j) - beta0(p) - mu(j) - beta1(p);
        nu_star = nu_eps(p) * n_rep(p) * w(p, j) + nu_alpha1(p);
        alpha_star = tmp * nu_eps(p) * n_rep(p) * w(p, j) /
                nu_star;
        if (gamma(p, i) == 1.0) {
            cur_val = alpha1(p, i);
            prop_val = 0.0;
            tmp = y_mean(p, j) - beta0(p) - mu(j) - beta1(p) -
                    cur_val;
            cur_den = //likelihood
                    -.5 * nu_eps(p) * w(p, j) * n_rep(p) * pow(tmp, 2) +
                    // prior for alpha1
                    log_omega +
                    .5 * log(nu_alpha1(p)) - .5 * (1.0 + nu_re) *
                    log1p(nu_alpha1(p) * pow(cur_val, 2) / nu_re);
            tmp = tmp + beta1(p) + cur_val - alpha0(p, i);
            prop_den = //likelihood
                    -.5 * nu_eps(p) * w(p, j) * n_rep(p) * pow(tmp, 2) +
                    // prior for alpha1
                    log1m_omega;

            cur_to_prop = 0.0;
            prop_to_cur = .5 * log(nu_star) - .5 * (1.0 + nu_re) *
                    log1p(nu_star * pow(cur_val - alpha_star, 2) / nu_re);
        } else {
            cur_val = 0.0;
            prop_val = RngStream_T(nu_re, rng) / sqrt(nu_star) +
                    alpha_star;

            tmp = y_mean(p, j) - beta0(p) - mu(j) - beta1(p) -
                    prop_val;
            prop_den = //likelihood
                    -.5 * nu_eps(p) * w(p, j) * n_rep(p) * pow(tmp, 2) +
                    // prior for alpha1
                    log_omega +
                    .5 * log(nu_alpha1(p)) - .5 * (1.0 + nu_re) *
                    log1p(nu_alpha1(p) * pow(prop_val, 2) / nu_re);
            tmp = tmp + beta1(p) + prop_val - alpha0(p, i);
            cur_den = //likelihood
                    -.5 * nu_eps(p) * w(p, j) * n_rep(p) * pow(tmp, 2) +
                    // prior for alpha1
                    log1m_omega;

            prop_to_cur = 0.0;
            cur_to_prop = .5 * log(nu_star) - .5 * (1.0 + nu_re) *
                    log1p(nu_star * pow(prop_val - alpha_star, 2) / nu_re);
        }
        // check proposal acceptance.
        if (update_alpha)
            alpha1(p, i) = log(RngStream_RandU01(rng)) <= prop_den + prop_to_cur -
                cur_den - cur_to_prop ? prop_val : cur_val;
        // "within-model" moves to re-sample weights, alpha1, alpha0
        if (alpha1(p, i) != 0.0) {
            if (update_gamma)
                gamma(p, i) = 1.0;
            j = 2 * i + 1;
            lambda_star = .5 * nu_re +
                    .5 * nu_alpha1(p) * pow(alpha1(p, i), 2);
            u(p, j) = RngStream_GA1(.5 + .5 * nu_re, rng) / lambda_star;
            nu_star = nu_eps(p) * n_rep(p) * w(p, j) +
                    nu_alpha1(p) * u(p, j);
            alpha_star = nu_eps(p) * n_rep(p) * w(p, j) *
                    (y_mean(p, j) - beta0(p) - mu(j) - beta1(p)) / nu_star;
            if (update_alpha)
                alpha1(p, i) = RngStream_N01(rng) / sqrt(nu_star) + alpha_star;

            j = 2 * i;
            nu_star = nu_eps(p) * n_rep(p) * w(p, j) +
                    nu_alpha0(p) * u(p, j);
            alpha_star = nu_eps(p) * n_rep(p) * w(p, j) *
                    (y_mean(p, j) - beta0(p) - mu(j)) / nu_star;
            if (update_alpha)
                alpha0(p, i) = RngStream_N01(rng) / sqrt(nu_star) + alpha_star;
        } else {
            // only sample alpha0, as alpha1 and u(p, 2 * i + 1) aren't
            // in model.
            if (update_gamma)
                gamma(p, i) = 0.0;
            nu_star = nu_eps(p) * n_rep(p) * (w(p, 2 * i) + w(p, 2 * i + 1)) +
                    nu_alpha0(p) * u(p, 2 * i);
            alpha_star = nu_eps(p) * n_rep(p) * (
                    w(p, 2 * i) * (y_mean(p, 2 * i) - beta0(p) - mu(2 * i)) +
                            w(p, 2*i + 1) * (y_mean(p, 2 * i + 1) - beta0(p) -
                                    mu(2 * i + 1))) / nu_star;
            if (update_alpha)
                alpha0(p, i) = RngStream_N01(rng) / sqrt(nu_star) + alpha_star;
        }
    }
}

void PairedMcmc::computeMahalaDist(const int & i) {
    double d;
    int j;
    const double mu_0(2 * i), mu_1(2 * i + 1);
    j = 2 * i
    for (int p = 0; p < n_peptide; p++) {
        d = nu_eps(p) * (y_ss(p, j) + n_rep(p) *
                pow(y_mean(p, j) - beta0(p) - mu_0 - alpha0(p, i), 2));
        mahala_dist(p, j) = d;
    }

    j = 2 * i + 1;
    for (int p = 0; p < n_peptide; p++) {
        if (gamma(p, i) == 0.0) {
            d = nu_eps(p) * (y_ss(p, j) + n_rep(p) *
                    pow(y_mean(p, j) - beta0(p) - mu_1 - alpha0(p, i), 2));
        } else {
            d = nu_eps(p) * (y_ss(p, j) + n_rep(p) *
                    pow(y_mean(p, j) - beta0(p) -
                            beta1(p) - mu_1 - alpha1(p, i), 2));
        }
        mahala_dist(p, j) = d;
    }
    return;
}

void PairedMcmc::updateWeights(const int & i, RngStream & rng) {
    // update the weights for alpha0, eps0, and eps1
    int p, j;
    double lambda_star, shape_star;

    j = 2 * i;
    for (p = 0; p < n_peptide; p++) {
        shape_star = .5 + .5 * nu_re;
        lambda_star = .5 * pow(alpha0(p, i), 2) * nu_alpha0(p) +
                .5 * nu_re;
        u(p, j) = RngStream_GA1(shape_star, rng) / lambda_star;

        shape_star = .5 * n_rep(p) + .5 * nu_err;
        lambda_star = .5 * mahala_dist(p, j) +
                .5 * (nu_err - 2.0);
        w(p, j) = RngStream_GA1(shape_star, rng) / lambda_star;
    }

    j = 2 * i + 1;
    for (p = 0; p < n_peptide; p++) {
        shape_star = .5 * n_rep(p) + .5 * nu_err;
        lambda_star = .5 * mahala_dist(p, j) +
                .5 * (nu_err - 2.0);
        w(p, j) = RngStream_GA1(shape_star, rng) / lambda_star;
    }
    return;
}

double PairedMcmc::dfJointDensity(const double x) {
    double lik(0.0);

    return(lik);
}






