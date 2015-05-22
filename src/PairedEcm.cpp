/*
 * PairedEcm.cpp
 *
 *  Created on: Nov 6, 2014
 *      Author: gimholte
 */

#include "RcppArmadillo.h"
#include <omp.h>
#include "MaximizeHypers.h"
#include "RngStream.h"
#include "Rv.h"
#include "utils.h"
#include "FitBase.h"
#include "PairedEcm.h"
#include "aghQuad.h"

PairedEcm::PairedEcm(const Rcpp::List & par_list, const Rcpp::List & data_list,
        const Rcpp::List & iter_list) :
        EcmAlgorithm(iter_list),
        m_beta0_prior_mean(Rcpp::as<double>(par_list["m_beta0_prior_mean"])),
        m_beta0_prior_prec(Rcpp::as<double>(par_list["m_beta0_prior_prec"])),
        m_beta1_prior_mean(Rcpp::as<double>(par_list["m_beta1_prior_mean"])),
        m_beta1_prior_prec(Rcpp::as<double>(par_list["m_beta1_prior_prec"])),
        s_alpha_rate(Rcpp::as<double>(par_list["s_alpha_rate"])),
        s_eps_rate(Rcpp::as<double>(par_list["s_eps_rate"])),
        lambda_alpha_rate(Rcpp::as<double>(par_list["lambda_alpha_rate"])),
        lambda_eps_rate(Rcpp::as<double>(par_list["lambda_eps_rate"])),
        nu_beta0_shape(Rcpp::as<double>(par_list["nu_beta0_shape"])),
        nu_beta0_rate(Rcpp::as<double>(par_list["nu_beta0_rate"])),
        nu_beta1_shape(Rcpp::as<double>(par_list["nu_beta1_shape"])),
        nu_beta1_rate(Rcpp::as<double>(par_list["nu_beta1_rate"])),
        lambda_a(Rcpp::as<double>(par_list["lambda_a"])),
        lambda_b(Rcpp::as<double>(par_list["lambda_b"])),
        nu(Rcpp::as<double>(par_list["nu"])),
        n_isamp_max(Rcpp::as<int>(iter_list["n_isamp"])) {
    // consume init_pars
    y_mean = Rcpp::as<arma::mat>(data_list["y_mean"]);
    y_ss = Rcpp::as<arma::mat>(data_list["y_ss"]);
    n_rep = Rcpp::as<arma::vec>(data_list["n_rep"]);
    pos = Rcpp::as<arma::ivec>(data_list["pos"]);
    n_position = Rcpp::as<int>(data_list["n_position"]);
    pos_start = Rcpp::as<arma::ivec>(data_list["pos_start"]);
    pos_count = Rcpp::as<arma::ivec>(data_list["pos_count"]);

    n_slide = y_mean.n_cols;
    n_subject = n_slide / 2;
    n_peptide = y_mean.n_rows;

    n_isamp_schedule = arma::ivec(n_iter);
    if (Rcpp::as<bool>(iter_list["schedule_iterations"])) {
        setImportanceSamplingSchedule();
    } else {
        n_isamp_schedule.fill(n_isamp_max);
    }
    n_isamp = n_isamp_schedule(0);

    beta0 = Rcpp::as<arma::vec>(par_list["beta0"]);
    beta1 = Rcpp::as<arma::vec>(par_list["beta1"]);
    m_beta0 = Rcpp::as<double>(par_list["m_beta0"]);
    m_beta1 = Rcpp::as<double>(par_list["m_beta1"]);
    nu_beta0 = Rcpp::as<double>(par_list["nu_beta0"]);
    nu_beta1 = Rcpp::as<double>(par_list["nu_beta1"]);

    mu = Rcpp::as<arma::vec>(par_list["mu"]);
    nu_eps = Rcpp::as<arma::vec>(par_list["nu_eps"]);
    nu_alpha0 = Rcpp::as<arma::vec>(par_list["nu_alpha0"]);
    nu_alpha1 = Rcpp::as<arma::vec>(par_list["nu_alpha1"]);
    s_eps = Rcpp::as<double>(par_list["s_eps"]);
    s_alpha = Rcpp::as<double>(par_list["s_alpha"]);
    lambda_eps = Rcpp::as<double>(par_list["lambda_eps"]);
    lambda_alpha = Rcpp::as<double>(par_list["lambda_alpha"]);
    a = Rcpp::as<arma::vec>(par_list["a"]);
    b = Rcpp::as<arma::vec>(par_list["b"]);
    omega = Rcpp::as<arma::vec>(par_list["omega"]);

    update_precision = Rcpp::as<bool>(iter_list["update_precision"]);
    update_precision_hypers = Rcpp::as<bool>(iter_list["update_precision_hypers"]);
    update_mu = Rcpp::as<bool>(iter_list["update_mu"]);
    update_beta0 = Rcpp::as<bool>(iter_list["update_beta0"]);
    update_beta1 = Rcpp::as<bool>(iter_list["update_beta1"]);
    update_beta_hypers = Rcpp::as<bool>(iter_list["update_beta_hypers"]);
    update_omega = Rcpp::as<bool>(iter_list["update_omega"]);
    update_weights = Rcpp::as<bool>(iter_list["update_weights"]);
    share_nu_alpha = false;//Rcpp::as<bool>(iter_list["share_nu_alpha"]);

    tau = arma::mat(n_peptide, n_subject);
    tau.fill(.001);

    // mu helpers
    mu_mean = arma::vec(n_slide, arma::fill::zeros);
    mu_temp = arma::vec(n_slide - 1, arma::fill::zeros);
    mu_prec_inner_diag = arma::vec(n_slide, arma::fill::zeros);
    mu_prec_inner = arma::mat(n_slide, n_slide, arma::fill::zeros);
    Q = Rcpp::as<arma::mat>(par_list["Q"]);
    mu_omega = arma::mat(n_slide - 1, n_slide - 1);

    // trace allocation
    hypers_tr = arma::mat(8, n_iter);
    beta0_tr = arma::mat(n_peptide, n_iter);
    beta1_tr = arma::mat(n_peptide, n_iter);
    mu_tr = arma::mat(n_slide, n_iter);
    nu_eps_tr = arma::mat(n_peptide, n_iter);
    nu_alpha0_tr = arma::mat(n_peptide, n_iter);
    nu_alpha1_tr = arma::mat(n_peptide, n_iter);
    omega_tr = arma::mat(n_peptide, n_iter);

    // t distribution weights
    w0_g0_samp = arma::cube(n_isamp_max, n_subject, n_peptide);
    w0_g0_samp.fill(1.0);

    w0_g1_samp = arma::cube(n_isamp_max, n_subject, n_peptide);
    w0_g1_samp.fill(1.0);

    w1_g0_samp = arma::cube(n_isamp_max, n_subject, n_peptide);
    w1_g0_samp.fill(1.0);

    w1_g1_samp = arma::cube(n_isamp_max, n_subject, n_peptide);
    w1_g1_samp.fill(1.0);

    u0_g0_samp = arma::cube(n_isamp_max, n_subject, n_peptide);
    u0_g0_samp.fill(1.0);

    u0_g1_samp = arma::cube(n_isamp_max, n_subject, n_peptide);
    u0_g1_samp.fill(1.0);

    u1_g1_samp = arma::cube(n_isamp_max, n_subject, n_peptide);
    u1_g1_samp.fill(1.0);

    q0_g0_weights = arma::cube(n_isamp_max, n_subject, n_peptide);
    q0_g0_weights.fill(1.0 / n_isamp_max);

    q0_g1_weights = arma::cube(n_isamp_max, n_subject, n_peptide);
    q0_g1_weights.fill(1.0 / n_isamp_max);

    q1_g1_weights = arma::cube(n_isamp_max, n_subject, n_peptide);
    q1_g1_weights.fill(1.0 / n_isamp_max);

    shrink_w0_g0_avg = arma::mat(n_peptide, n_subject);
    shrink_w1_g0_avg = arma::mat(n_peptide, n_subject);
    shrink_w0_g1_avg = arma::mat(n_peptide, n_subject);
    shrink_w1_g1_avg = arma::mat(n_peptide, n_subject);
    shrink_cross_g0_avg = arma::mat(n_peptide, n_subject);

    rng = std::vector<RngStream>(n_threads);
    gauss_work_x = std::vector<arma::vec>(n_threads);
    gauss_work_gx = std::vector<arma::vec>(n_threads);
    alpha_is_work = std::vector<arma::vec>(n_threads);
    density_inst_work = std::vector<arma::vec>(n_threads);
    density_joint_work = std::vector<arma::vec>(n_threads);
    gauss_rule_x = Rcpp::as<arma::vec>(data_list["gauss_rule_x"]);
    gauss_rule_w = Rcpp::as<arma::vec>(data_list["gauss_rule_w"]);
    n_rule = (int) gauss_rule_w.n_elem;

    for(int k = 0; k < n_threads; k++) {
        rng[k] = RngStream_CreateStream("");
        gauss_work_x[k].resize(n_rule);
        gauss_work_gx[k].resize(n_rule);
        alpha_is_work[k].resize(n_isamp_max);
        density_inst_work[k].resize(n_isamp_max);
        density_joint_work[k].resize(n_isamp_max);
    }
}

bool PairedEcm::checkInput() {
    bool is_ok = true;
    unsigned n_pep_us = (unsigned) n_peptide;
    unsigned n_sub_us = (unsigned) n_subject;
    is_ok &= y_ss.n_rows == y_mean.n_rows;
    is_ok &= y_ss.n_cols == y_mean.n_cols;
    is_ok &= n_pep_us == y_mean.n_rows;
    is_ok &= n_sub_us == (y_mean.n_cols / 2);
    is_ok &= (2 * n_sub_us) == y_mean.n_cols;
    is_ok &= nu_eps.n_elem == n_pep_us;
    is_ok &= nu_alpha0.n_elem == n_pep_us;
    is_ok &= beta0.n_elem == n_pep_us;
    is_ok &= beta1.n_elem == n_pep_us;
    is_ok &= mu.n_elem == (2 * n_sub_us);
    is_ok &= tau.n_rows == n_pep_us;
    is_ok &= tau.n_cols == n_sub_us;
    return is_ok;
}

void PairedEcm::setImportanceSamplingSchedule() {
    int n_break = n_iter > 10 ? 10 : n_iter;
    int n_per_break = (int) floor((double) n_iter / n_break);
    int n_isamp_min = n_iter > 5 ? 5 : n_iter;
    int b(1), v(n_isamp_min), j(0);
    // first break set has 1 iteration
    for (int i = 0; i < n_iter; i++) {
        n_isamp_schedule(i) = v;
        j++;
        if ((j == n_per_break) & (b < n_break)) {
            // reached a new break point, not the last one
            j = 0;
            b++;
            // ceil ok because b < n_break
            v = (int) ceil((n_isamp_max - n_isamp_min) *
                    (pow2((double) b / n_break))) + n_isamp_min;
        } else if ((j == n_per_break) & (b == n_break)) {
            v = n_isamp_max;
        }
    }
}

void PairedEcm::iterate(const int & iter_idx) {
    int p, i, th_id;
    n_isamp = n_isamp_schedule(iter_idx);
    #pragma omp parallel private(th_id, p) num_threads(n_threads)
    {
        th_id = omp_get_thread_num();
        #pragma omp for
        for (p = 0; p < n_peptide; p++) {
            if (update_weights) {
                eStep(p, th_id);
            }
            updateWeightSummaries(p);
        }
    }  // end omp parallel

    for (p = 0; p < n_peptide; p++) {
        if (update_omega)
            updateOmega(p);
    }

    for (p = 0; p < n_peptide; p++) {
        if (update_beta0)
            updateBeta0(p);
    }

    for (p = 0; p < n_peptide; p++) {
        if (update_beta1)
            updateBeta1(p);
    }

    for (i = 0; i < n_subject; i++) {
        prepMu(i);
    }

    if (update_mu)
        updateMu();
    if (update_beta_hypers)
        updateM0();
    if (update_beta_hypers)
        updateNuBeta0();
    if (update_beta_hypers)
        updateBeta1Hypers();
//    if (update_precision_hypers)
//        updatePrecHypers();

    for (p = 0; p < n_peptide; p++) {
        if (update_precision)
            updatePrecision(p);
    }
}

void PairedEcm::updateBeta0(const int & p) {
    double c_mean(0.0), c_prec(0.0);
    double d0, d1;
    int j;
    for (int i = 0; i < n_subject; i++) {
        j = 2 * i;
        d0 = y_mean(p, j) - mu(j);
        d1 = y_mean(p, j + 1) - mu(j + 1);
        c_mean += nu_eps(p) * n_rep(p) * nu_alpha0(p) * (1 - tau(p, i)) *
                (d0 * shrink_w0_g0_avg(p, i) +
                        d1 * shrink_w1_g0_avg(p, i)) ;
        c_mean += nu_eps(p) * n_rep(p) * nu_alpha1(p) * tau(p, i) *
                (d1 - beta1(p)) * shrink_w1_g1_avg(p, i);
        c_mean += nu_eps(p) * n_rep(p) * nu_alpha0(p) * tau(p, i) *
                d0 * shrink_w0_g1_avg(p, i);
        c_prec += nu_eps(p) * n_rep(p) * nu_alpha0(p) * (1 - tau(p, i)) *
                (shrink_w0_g0_avg(p, i) + shrink_w1_g0_avg(p, i));
        c_prec += nu_eps(p) * n_rep(p) * nu_alpha0(p) * tau(p, i) *
                shrink_w0_g1_avg(p, i) +
                nu_eps(p) * n_rep(p) * nu_alpha1(p) * tau(p, i) *
                shrink_w1_g1_avg(p, i);
    }
    c_mean += m_beta0 * nu_beta0;
    c_prec += nu_beta0;
    beta0(p) = c_mean / c_prec;
}

void PairedEcm::updateBeta1(const int & p) {
    double s(0.0), m(0.0), tmp;
    int j;
    for(int i = 0; i < n_subject; i++) {
            tmp = tau(p, i) * shrink_w1_g1_avg(p, i) *
                    n_rep(p) * nu_eps(p) * nu_alpha1(p);
            j = 2 * i + 1;
            m += tmp * (y_mean(p, j) - mu(j) - beta0(p));
            s += tmp;
    }

    m += m_beta1 * nu_beta1;
    s += nu_beta1;

    if (m < 0.0)
        beta1(p) = 0.0;
    else
        beta1(p) = m / s;
}

void PairedEcm::updateOmega(const int & p) {
    omega(p) = (a(pos(p)) - 1.0 + arma::sum(tau.row(p))) /
            (b(pos(p)) + a(pos(p)) - 2.0 + n_subject);
    if (omega(p) < 0.0) {
        omega(p) = 0.0;
    }
    if (omega(p) > 1.0) {
        omega(p) = 1.0;
    }
    return;
}

void PairedEcm::prepMu(const int & i) {
    double omega_i1_00(0.0), omega_i1_11(0.0);
    double omega_i0_00(0.0), omega_i0_11(0.0), omega_i0_01(0.0);
    double d_i1_0(0.0), d_i1_1(0.0), d_i0_0(0.0), d_i0_1(0.0);
    double m_prec, a_prec0, a_prec1, tmp;
    int j(2 * i);
    for (int p = 0; p < n_peptide; p++) {
        m_prec = n_rep(p) * nu_eps(p);
        a_prec0 = nu_alpha0(p);
        a_prec1 = nu_alpha1(p);

        omega_i0_00 += (1.0 - tau(p, i)) * m_prec *
                (a_prec0 * shrink_w0_g0_avg(p, i) +
                        m_prec * shrink_cross_g0_avg(p, i));
        omega_i0_11 += (1.0 - tau(p, i)) * m_prec *
                (a_prec0 * shrink_w1_g0_avg(p, i) +
                        m_prec * shrink_cross_g0_avg(p, i));
        omega_i0_01 += - (1.0 - tau(p, i)) * m_prec * m_prec *
                shrink_cross_g0_avg(p, i);

        omega_i1_00 += tau(p, i) * m_prec * a_prec0 *
                shrink_w0_g1_avg(p, i);
        omega_i1_11 += tau(p, i) * m_prec * a_prec1 *
                shrink_w1_g1_avg(p, i);

        tmp = (1.0 - tau(p, i));
        d_i0_0 += tmp * m_prec * a_prec0 * shrink_w0_g0_avg(p, i) *
                (y_mean(p, j) - beta0(p)) +
                tmp * m_prec * m_prec * shrink_cross_g0_avg(p, i) *
                (y_mean(p, j) - y_mean(p, j + 1));
        d_i0_1 += tmp * m_prec * a_prec0 * shrink_w1_g0_avg(p, i) *
                (y_mean(p, j + 1) - beta0(p)) +
                tmp * m_prec * m_prec * shrink_cross_g0_avg(p, i) *
                (y_mean(p, j + 1) - y_mean(p, j));

        tmp = tau(p, i) * m_prec;
        d_i1_0 += tmp * a_prec0 * shrink_w0_g1_avg(p, i) *
                (y_mean(p, j) - beta0(p));
        d_i1_1 += tmp * a_prec1 * shrink_w1_g1_avg(p, i) *
                (y_mean(p, j + 1) - beta0(p) - beta1(p));
    }

    mu_mean(j) = d_i0_0 + d_i1_0;
    mu_mean(j + 1) = d_i0_1 + d_i1_1;

    mu_prec_inner(j, j) = omega_i0_00 + omega_i1_00;
    mu_prec_inner(j, j + 1) = omega_i0_01;
    mu_prec_inner(j + 1, j) = omega_i0_01;
    mu_prec_inner(j + 1, j + 1) = omega_i0_11 + omega_i1_11;
}

void PairedEcm::updateMu() {
    // Get cholesky decomposition of conditional precision
    mu_omega = Q.t() * mu_prec_inner * Q;
    mu_omega.diag() += 1.0;
    mu_omega = arma::chol(mu_omega);
    // do forward/back solve to get the conditional mean
    mu_temp = arma::solve(mu_omega.t(), Q.t() * mu_mean);
    mu = Q * arma::solve(mu_omega, mu_temp);
}

void PairedEcm::updatePrecision(const int & p) {
    precPtrs prec_ptrs;
    prec_ptrs.y_mean = y_mean.memptr();
    prec_ptrs.y_ss = y_ss.memptr();
    prec_ptrs.mu = mu.memptr();
    prec_ptrs.tau = tau.memptr();

    prec_ptrs.q0_g0_weights = q0_g0_weights.slice(p).memptr();
    prec_ptrs.q0_g1_weights = q0_g1_weights.slice(p).memptr();
    prec_ptrs.q1_g1_weights = q1_g1_weights.slice(p).memptr();

    prec_ptrs.w0_g0_samp = w0_g0_samp.slice(p).memptr();
    prec_ptrs.w0_g1_samp = w0_g1_samp.slice(p).memptr();
    prec_ptrs.w1_g0_samp = w1_g0_samp.slice(p).memptr();
    prec_ptrs.w1_g1_samp = w1_g1_samp.slice(p).memptr();

    prec_ptrs.u0_g0_samp = u0_g0_samp.slice(p).memptr();
    prec_ptrs.u0_g1_samp = u0_g1_samp.slice(p).memptr();
    prec_ptrs.u1_g1_samp = u1_g1_samp.slice(p).memptr();

    prec_ptrs.n_subject = (int) n_subject;
    prec_ptrs.n_isamp = (int) n_isamp;
    prec_ptrs.n_isamp_max = (int) n_isamp_max;
    prec_ptrs.n_peptide = n_peptide;
    prec_ptrs.p = p;
    prec_ptrs.n_rep = n_rep(p);

    prec_ptrs.beta0 = beta0(p);
    prec_ptrs.beta1 = beta1(p);
    prec_ptrs.s_a = s_alpha;
    prec_ptrs.l_a = lambda_alpha;
    prec_ptrs.s_e = s_eps;
    prec_ptrs.l_e = lambda_eps;
    prec_ptrs.n_threads = n_threads;
    prec_ptrs.share_nu_alpha = share_nu_alpha;

    double prev_val[3] = {nu_eps(p), nu_alpha0(p), nu_alpha1(p)};
    double out[3];

    void * inputs = (void *) & prec_ptrs;
    maximizePrecision(inputs, prev_val, out, share_nu_alpha);

    nu_eps(p) = out[0];
    nu_alpha0(p) = out[1];
    nu_alpha1(p) = out[2];
}

void PairedEcm::updateM0() {
    m_beta0 = (m_beta0_prior_mean * m_beta0_prior_prec +
            arma::sum(beta0) * nu_beta0) /
                    (n_peptide * nu_beta0 + m_beta0_prior_prec);
}

void PairedEcm::updateNuBeta0() {
    double r(0.0), tmp;
    const double s0(nu_beta0_shape), l0(nu_beta0_rate);
    for(int p = 0; p < n_peptide; p++) {
        tmp = beta0(p) - m_beta0;
        r += pow(tmp, 2);
    }
    nu_beta0 = (n_peptide / 2.0 + (s0 - 1.0)) / (r / 2.0 + l0);
}

void PairedEcm::updatePrecHypers() {
    double par[2];
    if (share_nu_alpha) {
        maximizeGammaHypers(nu_alpha0.memptr(), nu_alpha0.n_elem,
                s_alpha_rate, lambda_alpha_rate, par);
    } else {
        arma::vec all_prec(n_peptide * 2);
        for (int r = 0; r < n_peptide; r++) {
            all_prec(2 * r) = nu_alpha0(r);
            all_prec(2 * r + 1) = nu_alpha1(r);
        }
        maximizeGammaHypers(all_prec.memptr(), all_prec.n_elem,
                s_alpha_rate, lambda_alpha_rate, par);
    }
    s_alpha = par[0];
    lambda_alpha = par[1];

    maximizeGammaHypers(nu_eps.memptr(), nu_eps.n_elem,
            s_eps_rate, lambda_eps_rate, par);
    s_eps = par[0];
    lambda_eps = par[1];

    return;
}

void PairedEcm::updateBeta1Hypers() {
    double par[2];
    par[0] = m_beta1;
    par[1] = nu_beta1;

    maximizeTnormHypers(beta1.memptr(), beta1.n_elem, m_beta1_prior_mean,
            m_beta1_prior_prec, /*nu_rate =*/ nu_beta1_rate,
            /*nu_shape =*/ nu_beta1_shape, par);

    m_beta1 = par[0];
    nu_beta1 = par[1];

    return;
}

void PairedEcm::eStep(const int & p, const int & th_id) {
    double den0;
    double den1[3];
    double offset0, offset1, delta0, delta1_g0, delta1_g1;
    double m0, m1, s0, s1, ap0, ap1;
    const double nu_e(nu_eps(p)), nu_a0(nu_alpha0(p)), nu_a1(nu_alpha1(p));
    const double b0(beta0(p)), b1(beta1(p));
    const double om(omega(p));

    int j0, j1;
    for (int i = 0; i < n_subject; i++) {
        j0 = 2 * i;
        j1 = j0 + 1;
        offset0 = y_ss(p, j0);
        offset1 = y_ss(p, j1);
        delta0 = y_mean(p, j0) - b0 - mu[j0];
        delta1_g0 = y_mean(p, j1) - b0 - mu[j1];
        delta1_g1 = delta1_g0 - b1;

        ap0 = n_rep(p) * nu_e + nu_a0;
        ap1 = n_rep(p) * nu_e + nu_a1;
        m0 = delta0 * n_rep(p) * nu_e / ap0;
        m1 = delta1_g1 * n_rep(p) * nu_e / ap1;
        s0 = 1.0 / sqrt(ap0);
        s1 = 1.0 / sqrt(ap1);

        den0 = densityGamma0(offset0, delta0, offset1, delta1_g0,
                nu, nu_e, nu_a0, n_rep[p],
                gauss_rule_x, gauss_rule_w,
                gauss_work_x[th_id], gauss_work_gx[th_id]);
        densityGamma1(den1, offset0, delta0, offset1, delta1_g1,
                nu, nu_e, nu_a0, nu_a1, n_rep[p],
                m0, m1, s0, s1,
                gauss_rule_x, gauss_rule_w,
                gauss_work_x[th_id], gauss_work_gx[th_id]);
        tau(p, i) = om * den1[2] / ((1 - om) * den0 + om * den1[2]);

        sampleWeightsGamma0(p, i, th_id, delta0, delta1_g0, offset0, offset1,
                alpha_is_work[th_id], density_inst_work[th_id],
                density_joint_work[th_id]);
        sampleWeightsGamma1(p, i, th_id, delta0, delta1_g1, offset0, offset1,
                alpha_is_work[th_id], density_inst_work[th_id],
                density_joint_work[th_id]);
    }
    return;
}

void PairedEcm::sampleWeightsGamma0(const int & p,
        const int & i, const int & th_id,
        const double & delta0, const double & delta1_g0,
        const double & offset0, const double & offset1,
        arma::vec & alpha_is_work, arma::vec & density_inst_work,
                    arma::vec & density_joint_work) {
    const double nu_e(nu_eps(p)), nu_a(nu_alpha0(p)), np(n_rep(p));
    double t, rate, shape, weightsum;

    // normalizing constant and scale/location of instrument density
    const double m = (delta0 + delta1_g0) * nu_e * np /
            (2 * nu_e * np + nu_a);
    const double sigma = sqrt(nu / ((nu - 2) * 2.0 * nu_e * np + nu_a));
    const double lg_const = R::lgammafn((1.0 + nu) / 2.0) -
            R::lgammafn(nu / 2.0) - .5 * log(nu * M_PI) - log(sigma);

    for (int r = 0; r < n_isamp; r++) {
        t = RngStream_T(nu, rng[th_id]);
        alpha_is_work[r] = t * sigma + m;
        density_inst_work[r] = lg_const - (.5 + nu / 2) *
                log1p(t * t / nu);
    }

    densityAlphaGamma0 g0(nu, np);
    g0.calc(density_joint_work, alpha_is_work,
            offset0, delta0, offset1, delta1_g0, nu, nu_e, nu_a, np);

    weightsum = 0.0;
    for (int r = 0; r < n_isamp; r++) {
        q0_g0_weights(r, i, p) = exp(density_joint_work[r] -
                density_inst_work[r]);
        weightsum += q0_g0_weights(r, i, p);

        shape = .5 * (nu + np);
        rate = .5 * nu + .5 *
                nu_e * (offset0 + np * pow(delta0 - alpha_is_work[r], 2));
        w0_g0_samp(r, i, p) = RngStream_GA1(shape, rng[th_id]) / rate;

        rate = .5 * nu + .5 *
                nu_e * (offset1 + np * pow(delta1_g0 - alpha_is_work[r], 2));
        w1_g0_samp(r, i, p) = RngStream_GA1(shape, rng[th_id]) / rate;

        shape = .5 * (nu + 1);
        rate = .5 * nu + .5 *
                nu_a * alpha_is_work[r] * alpha_is_work[r];
        u0_g0_samp(r, i, p) = RngStream_GA1(shape, rng[th_id]) / rate;
    }

    for (int r = 0; r < n_isamp; r++) {
        q0_g0_weights(r, i, p) = q0_g0_weights(r, i, p) / weightsum;
    }
    return;
}

void PairedEcm::sampleWeightsGamma1(const int & p,
        const int & i, const int & th_id,
        const double & delta0, const double & delta1_g1,
        const double & offset0, const double & offset1,
        arma::vec & alpha_is_work, arma::vec & density_inst_work,
        arma::vec & density_joint_work) {

    const double nu_e(nu_eps(p)), nu_a0(nu_alpha0(p)),
            nu_a1(nu_alpha1(p)), np(n_rep(p));
    double t, rate, weightsum;
    const double shape_1 = .5 * (1.0 + nu);
    const double shape_np = .5 * (np + nu);

    // normalizing constant and scale/location of instrument density
    const double m0 = delta0 * nu_e * np /
            (nu_e * np + nu_a0);
    const double m1 = delta1_g1 * nu_e * np /
            (nu_e * np + nu_a1);
    const double sigma0 = sqrt(nu / ((nu - 2) * 2.0 * nu_e * np + nu_a0));
    const double sigma1 = sqrt(nu / ((nu - 2) * 2.0 * nu_e * np + nu_a1));

    const double lg_const0 = R::lgammafn((1.0 + nu) / 2.0) -
            R::lgammafn(nu / 2.0) - .5 * log(nu * M_PI) - log(sigma0);
    const double lg_const1 = R::lgammafn((1.0 + nu) / 2.0) -
            R::lgammafn(nu / 2.0) - .5 * log(nu * M_PI) - log(sigma1);

    for (int r = 0; r < n_isamp; r++) {
        t = RngStream_T(nu, rng[th_id]);
        alpha_is_work[r] = t * sigma0 + m0;
        density_inst_work[r] = lg_const0 - (.5 + nu / 2) *
                log1p(t * t / nu);
    }

    densityAlphaGamma1 g1(nu, np);
    g1.calc(density_joint_work, alpha_is_work,
            offset0, delta0, nu, nu_e, nu_a0, np);
    weightsum = 0.0;
    for (int r = 0; r < n_isamp; r++) {
        q0_g1_weights(r, i, p) = exp(density_joint_work[r] -
                density_inst_work[r]);
        weightsum += q0_g1_weights(r, i, p);

        rate = .5 * nu + .5 *
                nu_e * (offset0 + np * pow(delta0 - alpha_is_work[r], 2));
        w0_g1_samp(r, i, p) = RngStream_GA1(shape_np, rng[th_id]) / rate;

        rate = .5 * nu + .5 *
                nu_a0 * alpha_is_work[r] * alpha_is_work[r];
        u0_g1_samp(r, i, p) = RngStream_GA1(shape_1, rng[th_id]) / rate;
    }
    for (int r = 0; r < n_isamp; r++) {
        q0_g1_weights(r, i, p) = q0_g1_weights(r, i, p) / weightsum;
    }

    for (int r = 0; r < n_isamp; r++) {
        t = RngStream_T(nu, rng[th_id]);
        alpha_is_work[r] = t * sigma1 + m1;
        density_inst_work[r] = lg_const1 - (.5 + nu / 2) *
                log1p(t * t / nu);
    }

    g1.calc(density_joint_work, alpha_is_work,
            offset1, delta1_g1, nu, nu_e, nu_a1, np);
    weightsum = 0.0;
    for (int r = 0; r < n_isamp; r++) {
        q1_g1_weights(r, i, p) = exp(density_joint_work[r] -
                density_inst_work[r]);
        weightsum += q1_g1_weights(r, i, p);

        rate = .5 * nu + .5 *
                nu_e * (offset1 + np * pow(delta1_g1 - alpha_is_work[r], 2));
        w1_g1_samp(r, i, p) = RngStream_GA1(shape_np, rng[th_id]) / rate;

        rate = .5 * nu + .5 *
                nu_a1 * alpha_is_work[r] * alpha_is_work[r];
        u1_g1_samp(r, i, p) = RngStream_GA1(shape_1, rng[th_id]) / rate;
    }
    for (int r = 0; r < n_isamp; r++) {
        q1_g1_weights(r, i, p) = q1_g1_weights(r, i, p) / weightsum;
    }
    return;
}

void PairedEcm::updateWeightSummaries(const int & p) {
    int i, r;
    double w0_g0, w1_g0, w0_g1, w1_g1;
    double u0_g0, u0_g1, u1_g1;
    double q, q0, q1;
    double bot, bot0, bot1;
    const double nu_e = nu_eps(p);
    const double nu_a0 = nu_alpha0(p);
    const double nu_a1 = nu_alpha1(p);

    double w0_g0_avg, w1_g0_avg, w0_g1_avg, w1_g1_avg, cross_avg;
    for (i = 0; i < n_subject; i++) {
        w0_g0_avg = 0.0; w0_g1_avg = 0.0;
        w1_g0_avg = 0.0; w1_g1_avg = 0.0;
        cross_avg = 0.0;
        for (r = 0; r < n_isamp; r++) {
            q = q0_g0_weights(r, i, p);
            q0 = q0_g1_weights(r, i, p);
            q1 = q1_g1_weights(r, i, p);
            w0_g0 = w0_g0_samp(r, i, p);
            w1_g0 = w1_g0_samp(r, i, p);
            w0_g1 = w0_g1_samp(r, i, p);
            w1_g1 = w1_g1_samp(r, i, p);
            u0_g0 = u0_g0_samp(r, i, p);
            u0_g1 = u0_g1_samp(r, i, p);
            u1_g1 = u1_g1_samp(r, i, p);

            bot = nu_e * n_rep(p) * (w0_g0 + w1_g0) +
                    nu_a0 * u0_g0;
            bot0 = nu_e * n_rep(p) * w0_g1 + nu_a0 * u0_g1;
            bot1 = nu_e * n_rep(p) * w1_g1 + nu_a1 * u1_g1;

            w0_g0_avg += q * w0_g0 * u0_g0 / bot;
            w1_g0_avg += q * w1_g0 * u0_g0 / bot;
            w0_g1_avg += q0 * w0_g1 * u0_g1 / bot0;
            w1_g1_avg += q1 * w1_g1 * u1_g1 / bot1;
            cross_avg += q * w0_g0 * w1_g0 / bot;
        }

        shrink_w0_g0_avg(p, i) = w0_g0_avg;
        shrink_w0_g1_avg(p, i) = w0_g1_avg;
        shrink_w1_g0_avg(p, i) = w1_g0_avg;
        shrink_w1_g1_avg(p, i) = w1_g1_avg;
        shrink_cross_g0_avg(p, i) = cross_avg;
    }
    return;
}

void PairedEcm::collectIteration(const int & iter_idx) {
    beta0_tr.col(iter_idx) = beta0;
    beta1_tr.col(iter_idx) = beta1;
    mu_tr.col(iter_idx) = mu;
    nu_eps_tr.col(iter_idx) = nu_eps;
    nu_alpha0_tr.col(iter_idx) = nu_alpha0;
    nu_alpha1_tr.col(iter_idx) = nu_alpha1;
    omega_tr.col(iter_idx) = omega;
    hypers_tr(0, iter_idx) = s_alpha;
    hypers_tr(1, iter_idx) = lambda_alpha;
    hypers_tr(2, iter_idx) = s_eps;
    hypers_tr(3, iter_idx) = lambda_eps;
    hypers_tr(4, iter_idx) = m_beta0;
    hypers_tr(5, iter_idx) = nu_beta0;
    hypers_tr(6, iter_idx) = m_beta1;
    hypers_tr(7, iter_idx) = nu_beta1;
    return;
}

const Rcpp::List PairedEcm::gatherOutput() {
    Rcpp::List out = Rcpp::List::create(
            Rcpp::Named("hypers") = hypers_tr,
            Rcpp::Named("beta0") = beta0_tr,
            Rcpp::Named("beta1") = beta1_tr,
            Rcpp::Named("mu") = mu_tr,
            Rcpp::Named("nu_eps") = nu_eps_tr,
            Rcpp::Named("nu_alpha0") = nu_alpha0_tr,
            Rcpp::Named("nu_alpha1") = nu_alpha1_tr,
            Rcpp::Named("omega") = omega_tr,
            Rcpp::Named("ppb") = tau);
    return(out);
}




