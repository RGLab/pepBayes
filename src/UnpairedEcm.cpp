#include "RcppArmadillo.h"
#include <omp.h>
#include "MaximizeHypers.h"
#include "RngStream.h"
#include "Rv.h"
#include "utils.h"
#include "FitBase.h"
#include "UnpairedEcm.h"
#include "aghQuad.h"

UnpairedEcm::UnpairedEcm(const Rcpp::List & par_list, const Rcpp::List & data_list,
        const Rcpp::List & iter_list) :
        EcmAlgorithm(iter_list),
        m_beta0_prior_mean(Rcpp::as<double>(par_list["m_beta0_prior_mean"])),
        m_beta0_prior_prec(Rcpp::as<double>(par_list["m_beta0_prior_prec"])),
        m_beta1_prior_mean(Rcpp::as<double>(par_list["m_beta1_prior_mean"])),
        m_beta1_prior_prec(Rcpp::as<double>(par_list["m_beta1_prior_prec"])),
        nu_beta0_shape(Rcpp::as<double>(par_list["nu_beta0_shape"])),
        nu_beta0_rate(Rcpp::as<double>(par_list["nu_beta0_rate"])),
        nu_beta1_shape(Rcpp::as<double>(par_list["nu_beta1_shape"])),
        nu_beta1_rate(Rcpp::as<double>(par_list["nu_beta1_rate"])),
        s_alpha_rate(Rcpp::as<double>(par_list["s_alpha_rate"])),
        s_eps_rate(Rcpp::as<double>(par_list["s_eps_rate"])),
        lambda_alpha_rate(Rcpp::as<double>(par_list["lambda_alpha_rate"])),
        lambda_eps_rate(Rcpp::as<double>(par_list["lambda_eps_rate"])),
        lambda_a(Rcpp::as<double>(par_list["lambda_a"])),
        lambda_b(Rcpp::as<double>(par_list["lambda_b"])),
        n_isamp_max(Rcpp::as<int>(iter_list["n_isamp"])),
        nu(Rcpp::as<double>(par_list["nu"]))

{
    // consume data_list
    y_mean_trt = Rcpp::as<arma::mat>(data_list["y_mean_trt"]);
    y_ss_trt = Rcpp::as<arma::mat>(data_list["y_ss_trt"]);
    y_mean_ctl = Rcpp::as<arma::mat>(data_list["y_mean_ctl"]);
    y_ss_ctl = Rcpp::as<arma::mat>(data_list["y_ss_ctl"]);

    n_rep = Rcpp::as<arma::vec>(data_list["n_rep"]);
    pos = Rcpp::as<arma::ivec>(data_list["pos"]);
    n_position = Rcpp::as<int>(data_list["n_position"]);
    pos_start = Rcpp::as<arma::ivec>(data_list["pos_start"]);
    pos_count = Rcpp::as<arma::ivec>(data_list["pos_count"]);

    n_peptide = y_mean_ctl.n_rows;
    n_trt = y_mean_trt.n_cols;
    n_ctl = y_mean_ctl.n_cols;
    n_slide = n_trt + n_ctl;

    n_isamp_schedule = arma::ivec(n_iter);
    if (Rcpp::as<bool>(iter_list["schedule_iterations"])) {
        setImportanceSamplingSchedule();
    } else {
        n_isamp_schedule.fill(n_isamp_max);
    }
    n_isamp = n_isamp_schedule(0);

    beta0 = Rcpp::as<arma::vec>(par_list["beta0"]);
    beta1 = Rcpp::as<arma::vec>(par_list["beta1"]);
    mu = Rcpp::as<arma::vec>(par_list["mu"]);

    m_beta0 = Rcpp::as<double>(par_list["m_beta0"]);
    m_beta1 = Rcpp::as<double>(par_list["m_beta1"]);
    nu_beta0 = Rcpp::as<double>(par_list["nu_beta0"]);
    nu_beta1 = Rcpp::as<double>(par_list["nu_beta1"]);

    nu_eps = Rcpp::as<arma::vec>(par_list["nu_eps"]);
    nu_alpha0 = Rcpp::as<arma::vec>(par_list["nu_alpha0"]);
    nu_alpha1 = Rcpp::as<arma::vec>(par_list["nu_alpha1"]);

    a = Rcpp::as<arma::vec>(par_list["a"]);
    b = Rcpp::as<arma::vec>(par_list["b"]);
    omega = Rcpp::as<arma::vec>(par_list["omega"]);

    s_eps = Rcpp::as<double>(par_list["s_eps"]);
    s_alpha = Rcpp::as<double>(par_list["s_alpha"]);
    lambda_eps = Rcpp::as<double>(par_list["lambda_eps"]);
    lambda_alpha = Rcpp::as<double>(par_list["lambda_alpha"]);
    tau = Rcpp::as<arma::mat>(par_list["gamma"]);

    w_ctl = arma::cube(n_isamp_max, n_ctl, n_peptide, arma::fill::ones);
    w_trt_g0 = arma::cube(n_isamp_max, n_trt, n_peptide, arma::fill::ones);
    w_trt_g1 = arma::cube(n_isamp_max, n_trt, n_peptide, arma::fill::ones);

    u_ctl = arma::cube(n_isamp_max, n_ctl, n_peptide, arma::fill::ones);
    u_trt_g0 = arma::cube(n_isamp_max, n_trt, n_peptide, arma::fill::ones);
    u_trt_g1 = arma::cube(n_isamp_max, n_trt, n_peptide, arma::fill::ones);

    q_ctl = arma::cube(n_isamp_max, n_ctl, n_peptide);
    q_ctl.fill(1.0 / n_isamp);
    q_trt_g0 = arma::cube(n_isamp_max, n_trt, n_peptide);
    q_trt_g0.fill(1.0 / n_isamp);
    q_trt_g1 = arma::cube(n_isamp_max, n_trt, n_peptide);
    q_trt_g1.fill(1.0 / n_isamp);

    shrink_wu_ctl_avg = arma::mat(n_peptide, n_ctl, arma::fill::ones);
    shrink_wu_g0_avg = arma::mat(n_peptide, n_trt, arma::fill::ones);
    shrink_wu_g1_avg = arma::mat(n_peptide, n_trt, arma::fill::ones);
    mean_shrink_g0 = arma::mat(n_peptide, n_trt, arma::fill::ones);
    mean_shrink_g1 = arma::mat(n_peptide, n_trt, arma::fill::ones);
    mean_shrink_ctl = arma::mat(n_peptide, n_ctl, arma::fill::ones);
    aprec_g0 = arma::mat(n_peptide, n_trt, arma::fill::ones);
    aprec_g1 = arma::mat(n_peptide, n_trt, arma::fill::ones);
    aprec_ctl = arma::mat(n_peptide, n_ctl, arma::fill::ones);

    // mu helpers
    mu_mean = arma::vec(n_slide, arma::fill::zeros);
    mu_temp = arma::vec(n_slide - 1, arma::fill::zeros);
    mu_prec_inner_diag = arma::vec(n_slide, arma::fill::zeros);
    mu_prec_inner = arma::mat(n_slide, n_slide, arma::fill::zeros);
    Q = Rcpp::as<arma::mat>(par_list["Q"]);
    mu_omega = arma::mat(n_slide - 1, n_slide - 1);

    update_precision = Rcpp::as<bool>(iter_list["update_precision"]);
    update_precision_hypers = Rcpp::as<bool>(iter_list["update_precision_hypers"]);
    update_mu = Rcpp::as<bool>(iter_list["update_mu"]);
    update_beta0 = Rcpp::as<bool>(iter_list["update_beta0"]);
    update_beta1 = Rcpp::as<bool>(iter_list["update_beta1"]);
    update_beta_hypers = Rcpp::as<bool>(iter_list["update_beta_hypers"]);
    update_omega = Rcpp::as<bool>(iter_list["update_omega"]);
    update_weights = Rcpp::as<bool>(iter_list["update_weights"]);
    share_nu_alpha = false;//Rcpp::as<bool>(iter_list["share_nu_alpha"]);

    // trace allocation
    hypers_tr = arma::mat(8, n_iter);
    beta0_tr = arma::mat(n_peptide, n_iter);
    beta1_tr = arma::mat(n_peptide, n_iter);
    mu_tr = arma::mat(n_slide, n_iter);
    nu_eps_tr = arma::mat(n_peptide, n_iter);
    nu_alpha0_tr = arma::mat(n_peptide, n_iter);
    nu_alpha1_tr = arma::mat(n_peptide, n_iter);
    omega_tr = arma::mat(n_peptide, n_iter);

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

void UnpairedEcm::setImportanceSamplingSchedule() {
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

void UnpairedEcm::updateWeightSummaries(const int & p) {
    const double nu_a0 = nu_alpha0(p);
    const double nu_a1 = nu_alpha1(p);
    const double nu_e = nu_eps(p);
    const double np = n_rep(p);

    double bot, bot0, bot1;
    double w, w0, w1;
    double u, u0, u1;
    double q, q0, q1;
    double shrink_ctl_avg, shrink_trt_g0_avg, shrink_trt_g1_avg;
    double mean_shrink_g0_tmp, mean_shrink_g1_tmp, mean_shrink_ctl_tmp;
    double aprec_g0_tmp, aprec_g1_tmp, aprec_ctl_tmp;
    int i, r;

    for (i = 0; i < n_ctl; i++) {
        shrink_ctl_avg = 0.0;
        mean_shrink_ctl_tmp = 0.0;
        aprec_ctl_tmp = 0.0;
        for (r = 0; r < n_isamp; r++) {
            w = w_ctl(r, i, p);
            u = u_ctl(r, i, p);
            q = q_ctl(r, i, p);
            bot = nu_e * n_rep(p) * w + nu_a0 * u;

            aprec_ctl_tmp += q * bot;
            shrink_ctl_avg += q * nu_e * np * w * u / bot;
            mean_shrink_ctl_tmp += q * nu_e * np * w / bot;
        }
        aprec_ctl(p, i) = aprec_ctl_tmp;
        shrink_wu_ctl_avg(p, i) = shrink_ctl_avg;
        mean_shrink_ctl(p, i) = mean_shrink_ctl_tmp;
    }

    for (i = 0; i < n_trt; i++) {
        shrink_trt_g0_avg = 0.0;
        shrink_trt_g1_avg = 0.0;
        mean_shrink_g0_tmp = 0.0;
        mean_shrink_g1_tmp = 0.0;
        aprec_g0_tmp = 0.0;
        aprec_g1_tmp = 0.0;
        for (r = 0; r < n_isamp; r++) {
            w0 = w_trt_g0(r, i, p);
            w1 = w_trt_g1(r, i, p);
            u0 = u_trt_g0(r, i, p);
            u1 = u_trt_g1(r, i, p);
            q0 = q_trt_g0(r, i, p);
            q1 = q_trt_g1(r, i, p);

            bot0 = nu_e * np * w0 + nu_a0 * u0;
            bot1 = nu_e * np * w1 + nu_a1 * u1;

            aprec_g0_tmp += q0 * bot0;
            aprec_g1_tmp += q1 * bot1;
            shrink_trt_g0_avg += q0 * nu_e * np * w0 * u0 / bot0;
            shrink_trt_g1_avg += q1 * nu_e * np * w1 * u1 / bot1;
            mean_shrink_g0_tmp += q0 * nu_e * np * w0 / bot0;
            mean_shrink_g1_tmp += q1 * nu_e * np * w1 / bot1;
        }

        aprec_g0(p, i) = aprec_g0_tmp;
        aprec_g1(p, i) = aprec_g1_tmp;
        shrink_wu_g0_avg(p, i) = shrink_trt_g0_avg;
        shrink_wu_g1_avg(p, i) = shrink_trt_g1_avg;
        mean_shrink_g0(p, i) = mean_shrink_g0_tmp;
        mean_shrink_g1(p, i) = mean_shrink_g1_tmp;
    }
}

bool UnpairedEcm::checkInput() {
    bool is_ok = true;
    unsigned n_pep_us = (unsigned) n_peptide;
    unsigned n_trt_us = (unsigned) n_trt;
    unsigned n_ctl_us = (unsigned) n_ctl;

    is_ok &= y_ss_trt.n_rows == y_mean_trt.n_rows;
    is_ok &= y_ss_trt.n_cols == y_mean_trt.n_cols;

    is_ok &= y_ss_ctl.n_rows == y_mean_ctl.n_rows;
    is_ok &= y_ss_ctl.n_cols == y_mean_ctl.n_cols;

    is_ok &= nu_eps.n_elem == n_pep_us;
    is_ok &= nu_alpha0.n_elem == n_pep_us;
    is_ok &= beta0.n_elem == n_pep_us;
    is_ok &= beta1.n_elem == n_pep_us;
    is_ok &= mu.n_elem == n_trt_us + n_ctl_us;

    return is_ok;
}

void UnpairedEcm::iterate(const int & iter_idx) {
    int p, i, th_id;
    n_isamp = n_isamp_schedule(iter_idx);
    #pragma omp parallel private(th_id, p) num_threads(n_threads)
    {
        th_id = omp_get_thread_num();
        #pragma omp for
        for (p = 0; p < n_peptide; p++) {
            eStep(p, th_id);
            updateWeightSummaries(p);
        }
    }

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

    for (i = 0; i < n_trt; i++) {
        prepMuTrt(i);
    }

    for (i = 0; i < n_ctl; i++) {
        prepMuCtl(i);
    }

    if (update_mu)
        updateMu();

    if (update_beta_hypers)
        updateM0();
    if (update_beta_hypers)
        updateNuBeta0();
    if (update_beta_hypers)
        updateBeta1Hypers();

    for (p = 0; p < n_peptide; p++) {
        if (update_precision)
            updatePrecision(p);
    }
//    if (update_precision_hypers)
//        updatePrecHypers();
}

void UnpairedEcm::updateBeta0(const int & p) {
    double c_mean(0.0), c_prec(0.0);
    double d, d0, d1;
    double s, s0, s1;
    for (int i = 0; i < n_ctl; i++) {
        s = nu_alpha0(p) * shrink_wu_ctl_avg(p, i);
        d = y_mean_ctl(p, i) - mu(i);
        c_mean += d * s;
        c_prec += s;
    }

    for (int i = 0; i < n_trt; i++) {
        s0 = nu_alpha0(p) * shrink_wu_g0_avg(p, i);
        s1 = nu_alpha1(p) * shrink_wu_g1_avg(p, i);

        d0 = y_mean_trt(p, i) - mu(n_ctl + i);
        d1 = d0 - beta1(p);

        c_mean += (1.0 - tau(p, i)) * s0 * d0;
        c_mean += tau(p, i) * s1 * d1;

        c_prec += (1.0 - tau(p, i)) * s0;
        c_prec += tau(p, i) * s1;
    }
    c_mean += m_beta0 * nu_beta0;
    c_prec += nu_beta0;

    beta0(p) = c_mean / c_prec;
}


void UnpairedEcm::updateBeta1(const int & p) {
    double s(0.0), m(0.0), tmp;
    int j;
    for(int i = 0; i < n_trt; i++) {
            tmp = tau(p, i) * nu_alpha1(p) * shrink_wu_g1_avg(p, i);
            j = n_ctl + i;
            m += tmp * (y_mean_trt(p, i) - mu(j) - beta0(p));
            s += tmp;
    }

    m += m_beta1 * nu_beta1;
    s += nu_beta1;
    if (m > 0) {
        beta1(p) =  m / s;
    } else {
        beta1(p) = 0.0;
    }
}

void UnpairedEcm::updateOmega(const int & p) {
    omega(p) = (a(pos(p)) - 1.0 + arma::sum(tau.row(p))) /
            (b(pos(p)) + a(pos(p)) - 2.0 + n_trt);
    if (omega(p) < 0.0) {
        omega(p) = 0.0;
    }
    if (omega(p) > 1.0) {
        omega(p) = 1.0;
    }
    return;
}

void UnpairedEcm::prepMuTrt(const int & i) {
    double omega_ii(0.0);
    double d(0.0);
    double tmp;
    for (int p = 0; p < n_peptide; p++) {
        omega_ii += (1 - tau(p, i)) *
                shrink_wu_g0_avg(p, i) * nu_alpha0(p);
        omega_ii += tau(p, i) *
                shrink_wu_g1_avg(p, i) * nu_alpha1(p);

        tmp = 1.0 - tau(p, i);
        d += tmp * shrink_wu_g0_avg(p, i) * nu_alpha0(p) *
                (y_mean_trt(p, i) - beta0(p));
        tmp = tau(p, i);
        d += tmp * shrink_wu_g1_avg(p, i) * nu_alpha1(p) *
                (y_mean_trt(p, i) - beta0(p) - beta1(p));
    }
    mu_mean(n_ctl + i) = d;
    mu_prec_inner(n_ctl + i, n_ctl + i) = omega_ii;
}

void UnpairedEcm::prepMuCtl(const int & i) {
    double omega_ii(0.0);
    double d(0.0);
    double a_prec;

    for (int p = 0; p < n_peptide; p++) {
        a_prec = nu_alpha0(p);
        omega_ii += shrink_wu_ctl_avg(p, i) * a_prec;
        d += shrink_wu_ctl_avg(p, i) * a_prec *
                (y_mean_ctl(p, i) - beta0(p));
    }

    mu_mean(i) = d;
    mu_prec_inner(i, i) = omega_ii;
}


void UnpairedEcm::updateMu() {
    // Get cholesky decomposition of conditional precision
    //mu_omega = Q.t() * arma::diagmat(mu_prec_inner_diag) * Q;
    mu_omega = Q.t() * mu_prec_inner * Q;
    mu_omega.diag() += 1.0;
    mu_omega = arma::chol(mu_omega);
    // do forward/back solve to get the conditional mean
    mu_temp = arma::solve(mu_omega.t(), Q.t() * mu_mean);
    mu = Q * arma::solve(mu_omega, mu_temp);
}

void UnpairedEcm::updatePrecision(const int & p) {
    precPtrsUnpaired prec_ptrs;
    prec_ptrs.y_mean_ctl = y_mean_ctl.memptr();
    prec_ptrs.y_mean_trt = y_mean_trt.memptr();
    prec_ptrs.y_ss_ctl = y_ss_ctl.memptr();
    prec_ptrs.y_ss_trt = y_ss_trt.memptr();
    prec_ptrs.mu = mu.memptr();
    prec_ptrs.tau = tau.memptr();

    prec_ptrs.w_ctl = w_ctl.slice(p).memptr();
    prec_ptrs.w_trt_g0 = w_trt_g0.slice(p).memptr();
    prec_ptrs.w_trt_g1 = w_trt_g1.slice(p).memptr();

    prec_ptrs.u_ctl = u_ctl.slice(p).memptr();
    prec_ptrs.u_trt_g0 = u_trt_g0.slice(p).memptr();
    prec_ptrs.u_trt_g1 = u_trt_g1.slice(p).memptr();

    prec_ptrs.q_ctl = q_ctl.slice(p).memptr();
    prec_ptrs.q_trt_g0 = q_trt_g0.slice(p).memptr();
    prec_ptrs.q_trt_g1 = q_trt_g1.slice(p).memptr();

    prec_ptrs.n_trt = (int) n_trt;
    prec_ptrs.n_ctl = (int) n_ctl;
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
    maximizePrecisionUnpaired(inputs, prev_val, out, share_nu_alpha);
    nu_eps(p) = out[0];
    nu_alpha0(p) = out[1];
    nu_alpha1(p) = out[2];
}

void UnpairedEcm::updateM0() {
    m_beta0 = (m_beta0_prior_mean * m_beta0_prior_prec +
            arma::sum(beta0) * nu_beta0) /
                    (n_peptide * nu_beta0 + m_beta0_prior_prec);
}

void UnpairedEcm::updateNuBeta0() {
    double r(0.0), tmp;
    const double s0(nu_beta0_shape), l0(nu_beta0_rate);
    for(int p = 0; p < n_peptide; p++) {
        tmp = beta0(p) - m_beta0;
        r += pow(tmp, 2);
    }
    nu_beta0 = (n_peptide / 2.0 + (s0 - 1.0)) / (r / 2.0 + l0);
}

void UnpairedEcm::updatePrecHypers() {
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

void UnpairedEcm::updateBeta1Hypers() {
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

void UnpairedEcm::eStep(const int & p, const int & th_id) {
    double den[3];
    double offset, delta, delta_g0, delta_g1;
    double m0, m1, s0, s1;
    const double nu_e(nu_eps(p)), nu_a0(nu_alpha0(p)), nu_a1(nu_alpha1(p));
    const double b0(beta0(p)), b1(beta1(p));
    const double om(omega(p));

    for(int i = 0; i < n_trt; i++){
        offset = y_ss_trt(p, i);
        delta_g0 = y_mean_trt(p, i) - b0 - mu(n_ctl + i);
        delta_g1 = delta_g0 - b1;
        m0 = delta_g0 * mean_shrink_g0(p, i);
        m1 = delta_g1 * mean_shrink_g1(p, i);
        s0 = 1 / sqrt(aprec_g0(p, i));
        s1 = 1.0 / sqrt(aprec_g1(p, i));

        densityGamma1(den, offset, delta_g0, offset, delta_g1,
                nu, nu_e, nu_a0, nu_a1, n_rep[p],
                m0, m1, s0, s1,
                gauss_rule_x, gauss_rule_w,
                gauss_work_x[th_id], gauss_work_gx[th_id]);

        if (update_weights) {
            tau(p, i) = om * den[1] / (om * den[1] + (1.0 - om) * den[0]);
        }
        if (update_weights) {
            sampleWeightsTrt(p, i, th_id, delta_g0, delta_g1, offset,
                    alpha_is_work[th_id],
                    density_inst_work[th_id],
                    density_joint_work[th_id]);
        }
    }

    for (int i = 0; i < n_ctl; i++) {
        offset = y_ss_ctl(p, i);
        delta = y_mean_ctl(p, i) - b0 - mu(i);
        if (update_weights) {
            sampleWeightsCtl(p, i, th_id, delta, offset,
                    alpha_is_work[th_id],
                    density_inst_work[th_id],
                    density_joint_work[th_id]);
        }
    }
    return;
}

void UnpairedEcm::sampleWeightsTrt(const int & p, const int & i,
        const int & th_id,
        const double & delta_g0, const double & delta_g1,
        const double & offset, arma::vec & alpha_is_work,
        arma::vec & density_inst_work,
        arma::vec & density_joint_work) {
    const double nu_e(nu_eps(p)), nu_a0(nu_alpha0(p)),
            nu_a1(nu_alpha1(p)), np(n_rep(p));
    double t, rate, weightsum;
    const double shape_1 = .5 * (1.0 + nu);
    const double shape_np = .5 * (np + nu);
    densityAlphaGamma1 g1(nu, np);

    // scale and location of instrument density
    const double m0 = delta_g0 * mean_shrink_g0(p, i);
    const double m1 = delta_g1 * mean_shrink_g1(p, i);
    const double sigma0 = 2.0 * sqrt(1 / aprec_g0(p, i));
    const double sigma1 = 2.0 * sqrt(1 / aprec_g1(p, i));

    const double df_samp(nu >= 2 ? nu - 1 : nu);
    const double lg_const0 = R::lgammafn((1.0 + df_samp) / 2.0) -
            R::lgammafn(df_samp) - .5 * log(df_samp * M_PI) - log(sigma0);
    const double lg_const1 = R::lgammafn((1.0 + nu) / 2.0) -
            R::lgammafn(df_samp / 2.0) - .5 * log(df_samp * M_PI) - log(sigma1);

    for (int r = 0; r < n_isamp; r++) {
        t = RngStream_T(df_samp, rng[th_id]);
        alpha_is_work[r] = t * sigma0 + m0;
        density_inst_work[r] = lg_const0 - (.5 + df_samp / 2) *
                log1p(t * t / df_samp);
    }

    g1.calc(density_joint_work, alpha_is_work,
            offset, delta_g0, nu, nu_e, nu_a0, np);
    weightsum = 0.0;
    for (int r = 0; r < n_isamp; r++) {
        q_trt_g0(r, i, p) = exp(density_joint_work[r] -
                density_inst_work[r]);
        weightsum += q_trt_g0(r, i, p);

        rate = .5 * nu + .5 *
                nu_e * (offset + np * pow(delta_g0 - alpha_is_work[r], 2));
        w_trt_g0(r, i, p) = RngStream_GA1(shape_np, rng[th_id]) / rate;

        rate = .5 * nu + .5 *
                nu_a0 * alpha_is_work[r] * alpha_is_work[r];
        u_trt_g0(r, i, p) = RngStream_GA1(shape_1, rng[th_id]) / rate;
    }
    // normalize weights
    for (int r = 0; r < n_isamp; r++) {
        q_trt_g0(r, i, p) = q_trt_g0(r, i, p) / weightsum;
    }

    for (int r = 0; r < n_isamp; r++) {
        t = RngStream_T(df_samp, rng[th_id]);
        alpha_is_work[r] = t * sigma1 + m1;
        density_inst_work[r] = lg_const1 - (.5 + df_samp / 2) *
                log1p(t * t / df_samp);
    }

    g1.calc(density_joint_work, alpha_is_work,
            offset, delta_g1, nu, nu_e, nu_a1, np);
    weightsum = 0.0;
    for (int r = 0; r < n_isamp; r++) {
        q_trt_g1(r, i, p) = exp(density_joint_work[r] -
                density_inst_work[r]);
        weightsum += q_trt_g1(r, i, p);

        rate = .5 * nu + .5 *
                nu_e * (offset + np * pow(delta_g1 - alpha_is_work[r], 2));
        w_trt_g1(r, i, p) = RngStream_GA1(shape_np, rng[th_id]) / rate;

        rate = .5 * nu + .5 *
                nu_a1 * alpha_is_work[r] * alpha_is_work[r];
        u_trt_g1(r, i, p) = RngStream_GA1(shape_1, rng[th_id]) / rate;
    }
    // normalize weights
    for (int r = 0; r < n_isamp; r++) {
        q_trt_g1(r, i, p) = q_trt_g1(r, i, p) / weightsum;
    }
    return;

}

void UnpairedEcm::sampleWeightsCtl(const int & p, const int & i, const int & th_id,
        const double & delta, const double & offset,
        arma::vec & alpha_is_work,
        arma::vec & density_inst_work,
        arma::vec & density_joint_work) {
    const double nu_e(nu_eps(p)), nu_a0(nu_alpha0(p)), np(n_rep(p));
     double t, rate, weightsum;
     const double shape_1 = .5 * (1.0 + nu);
     const double shape_np = .5 * (np + nu);
     densityAlphaGamma1 g1(nu, np);
     const double m0 = delta * mean_shrink_ctl(p, i);
     const double sigma0 = 2.0 * sqrt(1.0 / aprec_ctl(p, i));

     const double df_samp(nu >= 2 ? nu - 1 : nu);
     const double lg_const0 = R::lgammafn((1.0 + df_samp) / 2.0) -
             R::lgammafn(df_samp / 2.0) - .5 * log(df_samp * M_PI) - log(sigma0);

     for (int r = 0; r < n_isamp; r++) {
         t = RngStream_T(df_samp, rng[th_id]);
         alpha_is_work[r] = t * sigma0 + m0;
         density_inst_work[r] = lg_const0 - (.5 + df_samp / 2) *
                 log1p(t * t / df_samp);
     }

     g1.calc(density_joint_work, alpha_is_work,
             offset, delta, nu, nu_e, nu_a0, np);
     weightsum = 0.0;
     for (int r = 0; r < n_isamp; r++) {
         q_ctl(r, i, p) = exp(density_joint_work[r] -
                 density_inst_work[r]);
         weightsum += q_ctl(r, i, p);

         rate = .5 * nu + .5 *
                 nu_e * (offset + np * pow(delta - alpha_is_work[r], 2));
         w_ctl(r, i, p) = RngStream_GA1(shape_np, rng[th_id]) / rate;

         rate = .5 * nu + .5 *
                 nu_a0 * alpha_is_work[r] * alpha_is_work[r];
         u_ctl(r, i, p) = RngStream_GA1(shape_1, rng[th_id]) / rate;
     }
     // normalize weights
     for (int r = 0; r < n_isamp; r++) {
         q_ctl(r, i, p) = q_ctl(r, i, p) / weightsum;
     }

}

void UnpairedEcm::collectIteration(const int & iter_idx) {
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

const Rcpp::List UnpairedEcm::gatherOutput() {
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





