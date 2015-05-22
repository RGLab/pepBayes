/*
 * UnpairedEcm.cpp
 *
 *  Created on: Nov 21, 2014
 *      Author: gimholte
 */

#include "RcppArmadillo.h"
#include "RngStream.h"
#include "Rv.h"
#include <omp.h>
#include "Densities.h"
#include "AdaptiveGibbsSampler.h"
#include "FitBase.h"
#include "UnpairedMcmc.h"

UnpairedMcmc::UnpairedMcmc(const Rcpp::List & par_list,
        const Rcpp::List & data_list,
        const Rcpp::List & chain_list) : MarkovChain(chain_list),
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
        nu(Rcpp::as<double>(par_list["nu"])),
        m_beta1_tuner(n_burn),
        nu_beta1_tuner(n_burn)
{
    // set seed for arma's RNG generator so that initializations
    // are reproducible.
    std::srand(seed);

    y_mean_trt = Rcpp::as<arma::mat>(data_list["y_mean_trt"]);
    y_mean_control = Rcpp::as<arma::mat>(data_list["y_mean_ctl"]);
    y_ss_trt = Rcpp::as<arma::mat>(data_list["y_ss_trt"]);
    y_ss_control = Rcpp::as<arma::mat>(data_list["y_ss_ctl"]);

    n_rep = Rcpp::as<arma::ivec>(data_list["n_rep"]);
    pos = Rcpp::as<arma::ivec>(data_list["pos"]);
    pos_start = Rcpp::as<arma::ivec>(data_list["pos_start"]);
    pos_count = Rcpp::as<arma::ivec>(data_list["pos_count"]);

    n_position = pos_start.n_elem;
    n_peptide = y_ss_trt.n_rows;
    n_trt = y_ss_trt.n_cols;
    n_control = y_ss_control.n_cols;
    n_slide = n_trt + n_control;

    w_trt = Rcpp::as<arma::mat>(par_list["w_trt"]);
    w_ctl = Rcpp::as<arma::mat>(par_list["w_ctl"]);
    u_trt = Rcpp::as<arma::mat>(par_list["u_trt"]);
    u_ctl = Rcpp::as<arma::mat>(par_list["u_ctl"]);

    beta0 = Rcpp::as<arma::vec>(par_list["beta0"]);
    beta1 = Rcpp::as<arma::vec>(par_list["beta1"]);
    mu = Rcpp::as<arma::vec>(par_list["mu"]);
    alpha_trt = Rcpp::as<arma::mat>(par_list["alpha_trt"]);
    alpha_ctl = Rcpp::as<arma::mat>(par_list["alpha_ctl"]);

    gamma = Rcpp::as<arma::mat>(par_list["gamma"]);
    gamma_prob = arma::mat(n_peptide, n_trt);
    omega = Rcpp::as<arma::vec>(par_list["omega"]);
    a = Rcpp::as<arma::vec>(par_list["a"]);
    b = Rcpp::as<arma::vec>(par_list["b"]);

    nu_alpha0 = Rcpp::as<arma::vec>(par_list["nu_alpha0"]);
    nu_alpha1 = Rcpp::as<arma::vec>(par_list["nu_alpha1"]);
    nu_eps = Rcpp::as<arma::vec>(par_list["nu_eps"]);

    nu_beta0 = Rcpp::as<double>(par_list["nu_beta0"]);
    m_beta0 = Rcpp::as<double>(par_list["m_beta1"]);
    nu_beta1 = Rcpp::as<double>(par_list["nu_beta1"]);
    m_beta1 = Rcpp::as<double>(par_list["m_beta1"]);
    s_alpha = Rcpp::as<double>(par_list["s_alpha"]);
    lambda_alpha = Rcpp::as<double>(par_list["lambda_alpha"]);
    s_eps = Rcpp::as<double>(par_list["s_eps"]);
    lambda_eps = Rcpp::as<double>(par_list["lambda_eps"]);

    mu_mean = arma::vec(n_slide, arma::fill::zeros);
    mu_temp = arma::vec(n_slide - 1, arma::fill::zeros);
    mu_star_mean = arma::vec(n_slide - 1, arma::fill::zeros);
    mu_star = arma::vec(n_slide - 1, arma::fill::zeros);
    Q = Rcpp::as<arma::mat>(par_list["Q"]);
    mu_omega = arma::mat(n_slide - 1, n_slide - 1, arma::fill::zeros);
    mu_prec_inner_diag = arma::vec(n_slide, arma::fill::zeros);

    rng = std::vector<RngStream>(n_threads);
    for(int k = 0; k < n_threads; k++) {
        rng[k] = RngStream_CreateStream("");
    }

    s_alpha_sampler = AdaptiveGibbsSampler<ShapeConditional>(1.0, 20.0);
    s_eps_sampler = AdaptiveGibbsSampler<ShapeConditional>(1.0, 5.0);
    a_sampler = std::vector<AdaptiveGibbsSampler<BetaHyperConditional> >(n_position);
    b_sampler = std::vector<AdaptiveGibbsSampler<BetaHyperConditional> >(n_position);

    for(int q = 0; q < n_position; q++) {
        a_sampler[q] = AdaptiveGibbsSampler<BetaHyperConditional>(.1, 1.0);
        b_sampler[q] = AdaptiveGibbsSampler<BetaHyperConditional>(1.0, 10.0);
    }

    hypers_trace = arma::mat(8, n_samples);
    beta0_trace = arma::mat(n_peptide, n_samples);
    beta1_trace = arma::mat(n_peptide, n_samples);
    mu_trace = arma::mat(n_slide, n_samples);
    nu_alpha0_trace = arma::mat(n_peptide, n_samples);
    nu_alpha1_trace = arma::mat(n_peptide, n_samples);
    nu_eps_trace = arma::mat(n_peptide, n_samples);
    omega_trace = arma::mat(n_peptide, n_samples);
    a_trace = arma::mat(n_position, n_samples);
    b_trace = arma::mat(n_position, n_samples);
    ppa = arma::mat(n_peptide, n_trt, arma::fill::zeros);

    update_precision = Rcpp::as<bool>(chain_list["update_precision"]);
    update_precision_hypers = Rcpp::as<bool>(chain_list["update_precision_hypers"]);
    update_mu = Rcpp::as<bool>(chain_list["update_mu"]);
    update_beta0 = Rcpp::as<bool>(chain_list["update_beta0"]);
    update_beta1 = Rcpp::as<bool>(chain_list["update_beta1"]);
    update_beta_hypers = Rcpp::as<bool>(chain_list["update_beta_hypers"]);
    update_omega = Rcpp::as<bool>(chain_list["update_omega"]);
    update_omega_hypers = Rcpp::as<bool>(chain_list["update_omega_hypers"]);
    update_alpha = Rcpp::as<bool>(chain_list["update_alpha"]);
    update_weights = Rcpp::as<bool>(chain_list["update_weights"]);
    update_gamma = Rcpp::as<bool>(chain_list["update_gamma"]);

    Rcpp::Rcout << "Initialization complete. Running " <<
            (n_samples * n_thin + n_burn) << " total iterations" << std::endl;
    Rcpp::Rcout << "for " << n_slide << " slides with " << n_peptide <<
            " peptides." << std::endl;
}

void UnpairedMcmc::iterate() {
    int th_id, p, i, q;
    // the order in which updates occur is important
    #pragma omp parallel private(th_id, p, i, q) num_threads(n_threads)
    {
        th_id = omp_get_thread_num();
        #pragma omp for
        for (i = 0; i < n_trt; i++) {
            updateGammaAlphaWeights_Trt(i, rng[th_id]);
        }

        #pragma omp for
        for(i = 0; i < n_control; i++) {
            updateAlphaWeights_Ctl(i, rng[th_id]);
        }

        #pragma omp single
        {
            if (update_beta_hypers)
                updateM0(rng[th_id]);
        } // implied barrier

        #pragma omp for
        for (p = 0; p < n_peptide; p++) {
            if (update_beta1)
                updateBeta1(p, rng[th_id]);
        }

        #pragma omp for
        for (p = 0; p < n_peptide; p++) {
            if (update_beta0)
                updateBeta0(p, rng[th_id]);
        }

        #pragma omp for
        for (i = 0; i < n_trt; i++) {
            computeMuMean_Trt(i);
        }

        #pragma omp for
        for (i = 0; i < n_control; i++) {
            computeMuMean_Control(i);
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
                    updateNuBeta(rng[th_id]);
            }
        #pragma omp section
            {
                if (update_mu)
                    updateMu(rng[th_id]);
            }
        #pragma omp section
            {
                if (update_beta_hypers)
                    updateM1(rng[th_id]);
                if (update_beta_hypers)
                    updateNuBeta1(rng[th_id]);
            }
        } // end sections

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
            if (update_beta_hypers)
                updateBetaHypers(q, rng[th_id]);
        }

        #pragma omp for
        for (p = 0; p < n_peptide; p++) {
            if (update_omega)
                updateOmega(p, rng[th_id]);
        }
    }
}

void UnpairedMcmc::collectIteration(const int & sample_idx) {
    beta0_trace.col(sample_idx) = beta0;
    beta1_trace.col(sample_idx) = beta1;
    mu_trace.col(sample_idx) = mu;
    nu_alpha0_trace.col(sample_idx) = nu_alpha0;
    nu_alpha1_trace.col(sample_idx) = nu_alpha1;
    nu_eps_trace.col(sample_idx) = nu_eps;
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

    ppa = ppa * ((double) sample_idx / (sample_idx + 1.0)) + gamma_prob / (sample_idx + 1.0);
}

const Rcpp::List UnpairedMcmc::chainOutput() {
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
            // The remaining parameters are returned for chain restarting
            Rcpp::Named("u_trt") = u_trt,
            Rcpp::Named("u_ctl") = u_ctl,
            Rcpp::Named("w_trt") = w_trt,
            Rcpp::Named("w_ctl") = w_ctl,
            Rcpp::Named("gamma") = gamma,
            Rcpp::Named("alpha_trt") = alpha_trt,
            Rcpp::Named("alpha_ctl") = alpha_ctl);
    return out;
}

void UnpairedMcmc::updateGammaAlphaWeights_Trt(const int & i,
        RngStream & rng) {
    int p;
    const double mu1 = mu(n_control + i);

    double delta0, delta1;
    double shrink0, shrink1;
    double nu_star0, nu_star1;
    double mu_star0, mu_star1;
    double lp0, lp1, prob;

    double alpha_mean, alpha_prec, lambda_w, lambda_u, delta;

    for(p = 0; p < n_peptide; p++) {
        delta0 = (y_mean_trt(p, i) - beta0(p) - mu1);
        delta1 = delta0 - beta1(p);

        nu_star0 = nu_eps(p) * n_rep(p) * w_trt(p, i) +
                nu_alpha0(p) * u_trt(p, i);
        nu_star1 = nu_eps(p) * n_rep(p) * w_trt(p, i) +
                nu_alpha1(p) * u_trt(p, i);
        shrink0 = nu_eps(p) * n_rep(p) * w_trt(p, i) / nu_star0;
        shrink1 = nu_eps(p) * n_rep(p) * w_trt(p, i) / nu_star1;

        mu_star0 = delta0 * shrink0;
        mu_star1 = delta1 * shrink1;

        lp0 = log1mFromLogit(omega(p));
        lp1 = logFromLogit(omega(p));

        lp0 += .5 * log(nu_alpha0(p)) - .5 * log(nu_star0);
        lp1 += .5 * log(nu_alpha1(p)) - .5 * log(nu_star1);

        lp0 += -.5 * u_trt(p, i) * nu_alpha0(p) * delta0 * delta0 * shrink0;
        lp1 += -.5 * u_trt(p, i) * nu_alpha1(p) * delta1 * delta1 * shrink1;

        prob = 1.0 / (1.0 + exp(lp0 - lp1));
        gamma_prob(p, i) = prob;
        if (update_gamma) {
            gamma(p, i) = (RngStream_RandU01(rng) <= prob) ? 1.0 : 0.0;
        }

        if ((gamma(p, i) > 0.0) ) {
            alpha_mean = mu_star1;
            alpha_prec = nu_star1;
        } else {
            alpha_mean = mu_star0;
            alpha_prec = nu_star0;
        }

        if (update_alpha) {
            alpha_trt(p, i) = RngStream_N01(rng) / sqrt(alpha_prec) + alpha_mean;
        }

        if (gamma(p, i) > 0.0) {
            delta = delta1 - alpha_trt(p, i);
        } else {
            delta = delta0 - alpha_trt(p, i);
        }

        if (update_weights) {
            lambda_w = .5 * nu + .5 * nu_eps(p) * n_rep(p) * delta * delta +
                    .5 * nu_eps(p) * y_ss_trt(p, i);
            w_trt(p, i) = RngStream_GA1(.5 * nu + .5 * n_rep(p), rng) / lambda_w;

            lambda_u = gamma(p, i) > 0.0 ? nu_alpha1(p) : nu_alpha0(p);
            lambda_u *= .5 * alpha_trt(p, i) * alpha_trt(p, i);
            lambda_u += .5 * nu;
            u_trt(p, i) = RngStream_GA1(.5 * nu, rng) / lambda_u;
        }
    }
    return;
}

void UnpairedMcmc::updateAlphaWeights_Ctl(const int & i, RngStream & rng) {
    int p;
    double mu0 = mu(i);
    double nu_star;
    double mu_star;
    double delta;
    double lambda;
    for(p = 0; p < n_peptide; p++) {
        delta = y_mean_control(p, i) - beta0(p) - mu0;
        nu_star = n_rep(p) * nu_eps(p) * w_ctl(p, i) +
                nu_alpha0(p) * u_ctl(p, i);
        mu_star = n_rep(p) * nu_eps(p) * w_ctl(p, i) * delta / nu_star;
        if (update_alpha) {
            alpha_ctl(p, i) = RngStream_N01(rng) / sqrt(nu_star) + mu_star;
        }
        delta = delta - alpha_ctl(p, i);
        lambda = .5 * nu_eps(p) *
                (n_rep(p) * delta * delta + y_ss_control(p, i)) +
                .5 * nu;
        if (update_weights) {
            w_ctl(p, i) = RngStream_GA1(.5 * (n_rep(p) + nu), rng) / lambda;
        }
        lambda = .5 * nu_alpha0(p) * alpha_ctl(p, i) * alpha_ctl(p, i) +
                .5 * nu;
        if (update_weights) {
            u_ctl(p, i) = RngStream_GA1(.5 * (1.0 + nu), rng) / lambda;
        }
    }
    return;
}

void UnpairedMcmc::computeMuMean_Trt(const int & i) {
    mu_mean(n_control + i) = arma::sum(
            (y_mean_trt.col(i) - beta0 - alpha_trt.col(i) -
                    gamma.col(i) % beta1) % nu_eps % n_rep % w_trt.col(i));

    mu_prec_inner_diag(n_control + i) = arma::sum(nu_eps % n_rep % w_trt.col(i));
}

void UnpairedMcmc::computeMuMean_Control(const int & i) {
    mu_mean(i) = arma::sum(
            (y_mean_control.col(i) - beta0 -
                    alpha_ctl.col(i)) % nu_eps % n_rep % w_ctl.col(i));

    mu_prec_inner_diag(i) = arma::sum(nu_eps % n_rep % w_ctl.col(i));
}

void UnpairedMcmc::updateMu(RngStream & rng) {
    // Get cholesky decomposition of conditional precision
    mu_omega = Q.t() * arma::diagmat(mu_prec_inner_diag) * Q;
    mu_omega.diag() += 1.0;
    mu_omega = arma::chol(mu_omega);
    // do forward/back solve to get the conditional mean
    mu_temp = arma::solve(mu_omega.t(), Q.t() * mu_mean);
    mu_star_mean = arma::solve(mu_omega, mu_temp);
    // fill with n(0, 1) values
    for(int i = 0; i < n_slide - 1; i++) {
        mu_temp(i) = RngStream_N01(rng);
    }
    // solve linear system to get N(0, omega^(-1)) multivariate.
    mu_star = arma::solve(mu_omega, mu_temp) + mu_star_mean;
    // transform mu_star to recover mu
    mu = Q * mu_star;
}

void UnpairedMcmc::updateM0(RngStream & rng) {
    double m(arma::sum(beta0)), s(0.0);
    m *= nu_beta0;
    m += m_beta0_prior_mean * m_beta0_prior_prec;
    s = nu_beta0 * n_peptide + m_beta0_prior_prec;
    m /= s;
    m_beta0 = RngStream_N01(rng) / sqrt(s) + m;
}

void UnpairedMcmc::updateBeta0(const int & p, RngStream & rng) {
    double s(0.0), m(0.0), u, prec(n_rep(p) * nu_eps(p));
    for(int i = 0; i < n_trt; i++) {
        u = w_trt(p, i);
        m += u * prec * (y_mean_trt(p, i) - alpha_trt(p, i) - mu(n_control + i) -
                gamma(p, i) * beta1(p));
        s += u * prec;
    }
    for(int j = 0; j < n_control; j++) {
        u = w_ctl(p, j);
        m += u * prec * (y_mean_control(p, j) - alpha_ctl(p, j) - mu(j));
        s += u * prec;
    }
    m += m_beta0 * nu_beta0;
    s += nu_beta0;
    beta0(p) = RngStream_N01(rng) / sqrt(s) + (m / s);
}

void UnpairedMcmc::updateBeta1(const int & p, RngStream & rng) {
    double s(0.0), m(0.0), tmp;
    for(int i = 0; i < n_trt; i++) {
        if (gamma(p, i) > 0.0) {
            tmp = w_trt(p, i);
            m += tmp * (y_mean_trt(p, i) - mu(n_control + i) - beta0(p) -
                    alpha_trt(p, i));
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

void UnpairedMcmc::updateAlphaHypers(RngStream & rng) {
    int p;
    double density_args[4];
    double log_sum(0.0);
    double nu_sum(0.0);
    for(p = 0; p < n_peptide; p++) {
        nu_sum += nu_alpha0(p);
        log_sum += log(nu_alpha0(p));
        nu_sum += nu_alpha1(p);
        log_sum += log(nu_alpha1(p));
    }
    density_args[0] = (double) 2.0 * n_peptide;
    density_args[1] = log(lambda_alpha_rate + nu_sum);
    density_args[2] = log_sum;
    density_args[3] = s_alpha_rate;
    if(!R_FINITE(log_sum)) {
        Rcpp::stop("invalid input to alpha hyper sampler");
    }
    s_alpha = s_alpha_sampler.sample(rng, density_args);
    lambda_alpha = RngStream_GA1(1.0 + 2.0 * n_peptide * s_alpha, rng) /
            (lambda_alpha_rate + nu_sum);
}

void UnpairedMcmc::updateEpsilonHypers(RngStream & rng) {
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
    density_args[3]  = s_eps_rate;
    s_eps = s_eps_sampler.sample(rng, density_args);
    lambda_eps = RngStream_GA1(1.0 + n_peptide * s_eps, rng) /
            (lambda_eps_rate + nu_sum);
}

void UnpairedMcmc::updateNuBeta(RngStream & rng) {
    int p;
    double s(0.0);
    for(p = 0; p < n_peptide; p++) {
        s += (beta0(p) - m_beta0) * (beta0(p) - m_beta0);
    }
    s /= 2.0;
    s += nu_beta0_rate;
    nu_beta0 = RngStream_GA1(nu_beta0_rate + n_peptide / 2.0, rng) / s;
}

void UnpairedMcmc::updateNuEps(const int & p, RngStream & rng) {
    double d0, s(0.0);
    const double np(n_rep(p));
    for(int i = 0; i < n_control; i++) {
        d0 = y_mean_control(p, i) - beta0(p) - mu(i) - alpha_ctl(p, i);
        d0 *= d0 * np;
        d0 += y_ss_control(p, i);
        d0 *= w_ctl(p, i);
        s += d0;
    }

    for(int i = 0; i < n_trt; i++) {
        d0 = y_mean_trt(p, i) - beta0(p) - mu(n_control + i) -
                alpha_trt(p, i) - gamma(p, i) * beta1(p);
        d0 *= d0 * np;
        d0 += y_ss_trt(p, i);
        d0 *= w_trt(p, i);
        s += d0;
    }

    s /= 2.0;
    s += lambda_eps;
    nu_eps(p) = RngStream_GA1(s_eps + n_slide * np / 2.0, rng) / s;
}

void UnpairedMcmc::updateNuAlpha(const int & p, RngStream & rng) {
    double s0(0.0), s1(0.0), a0, ell(0.0);
    int i;
    for(i = 0; i < n_control; i++) {
        a0 = alpha_ctl(p, i);
        s0 += a0 * a0 * u_ctl(p, i);
    }

    for(i = 0; i < n_trt; i++) {
        a0 = alpha_trt(p, i);
        if (gamma(p, i) > 0.0) {
            s1 += a0 * a0 * u_trt(p, i);
            ell += 1.0;
        } else {
            s0 += a0 * a0 * u_trt(p, i);
        }
    }
    s0 /= 2.0;
    s0 += lambda_alpha;
    s1 /= 2.0;
    s1 += lambda_alpha;
    nu_alpha0(p) = RngStream_GA1(s_alpha + (n_control + n_trt - ell) / 2.0,
            rng) / s0;
    nu_alpha1(p) = RngStream_GA1(s_alpha + ell / 2.0,
            rng) / s1;
}

void UnpairedMcmc::updateOmega(const int & p, RngStream & rng) {
    const double s = arma::sum(gamma.row(p));
    omega(p) = RngStream_LogitBeta(a(pos(p)) + s, b(pos(p)) + n_trt - s, rng);
}

void UnpairedMcmc::updateM1(RngStream & rng) {
    const double sum_b1 = arma::sum(beta1);
    const double cond_tau = n_peptide * nu_beta1 + m_beta1_prior_prec;
    const double cond_mu =  nu_beta1 * sum_b1 +
            m_beta1_prior_mean * m_beta1_prior_prec;
    const double prop_m1 = m_beta1 * exp((RngStream_RandU01(rng) - .5) *
            m_beta1_tuner.getScale());

    const double cur_lik = - m_beta1 * m_beta1 * cond_tau / 2 +
            cond_mu * m_beta1 -
            n_peptide * R::pnorm(0.0, m_beta1, 1.0 / sqrt(nu_beta1), 0, 1);
    const double prop_lik = - prop_m1 * prop_m1 * cond_tau / 2 +
            cond_mu * prop_m1 -
            n_peptide * R::pnorm(0.0, prop_m1, 1.0 / sqrt(nu_beta1), 0, 1);

    const double cur_to_prop = - log(prop_m1);
    const double prop_to_cur = - log(m_beta1);

    m_beta1 = metropolisHastings(m_beta1, cur_lik, cur_to_prop,
            prop_m1, prop_lik, prop_to_cur, rng);
    m_beta1_tuner.update(m_beta1 == prop_m1);
    return;
}

void UnpairedMcmc::updateNuBeta1(RngStream & rng) {
    const double ss_beta1 = pow(arma::norm(beta1 - m_beta1, 2), 2);
    const double cond_shape = n_peptide / 2.0 + .5;
    const double cond_rate = ss_beta1 / 2.0 + .5;
    const double prop_nb1 = nu_beta1 * exp((RngStream_RandU01(rng) - .5) *
            nu_beta1_tuner.getScale());

    const double cur_lik = - nu_beta1 * cond_rate +
            cond_shape * log(nu_beta1) -
            n_peptide * R::pnorm(0.0, m_beta1,
                    1.0 / sqrt(nu_beta1), 0, 1);
    const double prop_lik = -prop_nb1 * cond_rate +
            (cond_shape - 1.0) * log(prop_nb1) -
            n_peptide * R::pnorm(0.0, m_beta1,
                    1.0 / sqrt(prop_nb1), 0, 1);

    const double cur_to_prop = -log(prop_nb1);
    const double prop_to_cur = -log(nu_beta1);

    nu_beta1 = metropolisHastings(nu_beta1, cur_lik, cur_to_prop, prop_nb1,
            prop_lik, prop_to_cur, rng);
    nu_beta1_tuner.update(nu_beta1 == prop_nb1);
}

void UnpairedMcmc::updateBetaHypers(const int & q, RngStream & rng) {
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







