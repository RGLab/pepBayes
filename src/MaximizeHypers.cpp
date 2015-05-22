/*
 * MaximizeHypers.cpp
 *
 *  Created on: Jun 27, 2014
 *      Author: gimholte
 */

#include "R.h"
#include "Rmath.h"
#include <R_ext/Applic.h>
#include "MaximizeHypers.h"
#include "RngStream.h"
#include "utils.h"
#include <cfloat>

void maximizeGammaHypers(const double * data, int data_len, const double l_s,
        const double l_lambda, double * x) {
    // inside gammaHyperObjectiveFn:
    // const double sum_log_x(input[0]);
    // const double sum_x(input[1]);
    // const double P(input[2]);
    // const double l_s(input[3]);
    // const double l_l(input[4]);

    double input[5];
    input[0] = 0.0;
    input[1] = 0.0;
    for(int j = 0; j < data_len; j++) {
        input[0] += log(data[j]);
        input[1] += data[j];
    }
    input[2] = (double) data_len;
    input[3] = l_s;
    input[4] = l_lambda;

    // initial guess for max of shape (pretty good, actually)
    double par, s = log(input[1] / input[2]) - input[0] / input[2];
    s = ((3 - s) + sqrt((s - 3) * (s - 3) + 24.0 * s)) / (12 * s);
    par = s;

    // misc BFGS tuning parameters and reporting parameters
    int n(1), mask(1), report(10), maxit(100), fncount, grcount, fail, trace(0);
    double Fmin, abstol(- INFINITY), reltol(1.0e-8);

    vmmin(n, &par, &Fmin, gammaHyperObjectiveFn, gammaHyperObjectiveGr,
            maxit, trace, &mask, abstol, reltol, report, (void *) input,
            &fncount, &grcount, &fail);

    // optimization failed for some reason; return good guess
    if (fail > 0)
        par = s;

    x[0] = par;
    x[1] = input[2] * par / (l_lambda + input[1]);
    return;
}

void maximizeTnormHypers(const double * data, int data_len, const double m_mean,
        const double m_prec, const double nu_rate, const double nu_shape,
        double * x) {
    // inside betaHyperObjectiveFn:
    // const double beta1_sum(input[0]);
    // const double beta1_sqr_sum(input[1]);
    // const double P(input[2]);
    // const double l1(input[3]);
    // const double s1(input[4]);
    // const double m_bar(input[5]);
    // const double nu_m(input[6]);

    double input[7];
    input[0] = 0.0;
    input[1] = 0.0;
    for (int j = 0; j < data_len; j++) {
        input[0] += data[j];
        input[1] += data[j] * data[j];
    }
    input[2] = (double) data_len;
    input[3] = nu_rate;
    input[4] = nu_shape;
    input[5] = m_mean;
    input[6] = m_prec;

    // initial guess for location and (log) inverse-scale;
    double par[2];
    par[0] = input[0] / input[2];
    par[1] = - log(1.0 + input[1] / input[2] - par[0] * par[0]);

    // misc BFGS tuning parameters and reporting parameters
    int n(2), report(10), maxit(100), fncount, grcount, fail, trace(0);
    int mask[2] = {1, 1};
    double Fmin, abstol(- INFINITY), reltol(1.0e-8);

    vmmin(n, par, &Fmin, betaHyperObjectiveFn, betaHyperObjectiveGr,
            maxit, trace, mask, abstol, reltol, report, (void *) input,
            &fncount, &grcount, &fail);

    if (fail > 0)
        return;

    x[0] = par[0];
    x[1] = exp(par[1]);
}


void maximizePrecision(void * inputs, double * prev_val, double * x,
        bool share_nu_alpha) {
    precPtrs * ex = static_cast<precPtrs * >(inputs);
    int p = ex->p;
    int n_pep = ex->n_peptide;
    int n_subject = ex->n_subject;
    double TAU_CHECK_VAL = 1.0e-15;
    double tau_sum(0.0);
    bool orig_share = share_nu_alpha;

    for (int i = 0; i < n_subject; i++) {
        tau_sum += ex->tau[n_pep * i + p];
    }
    if ((tau_sum < TAU_CHECK_VAL) && !share_nu_alpha) {
        ex->share_nu_alpha = true;
        share_nu_alpha = true;
    }

    // starting location
    double par[3];
    double par_out[3];
    int n;
    par[0] = log(prev_val[0]);
    par[1] = log(prev_val[1]);
    if (share_nu_alpha) {
        par[2] = log(prev_val[1]);
        n = 2;
    } else {
        par[2] = log(prev_val[2]);
        n = 3;
    }

    // misc BFGS tuning parameters and reporting parameters
    int maxit(100), fncount, fail, trace(0);
    double Fmin, abstol(- INFINITY), intol(1.0e-4);
    double alpha(1.0), beta(.5), gamma(2.0);

    nmmin(n, par, par_out, &Fmin, precisionObjectiveFn,
            &fail, abstol, intol, inputs, alpha, beta, gamma, trace,
            &fncount, maxit);

    if (fail > 0) {
        x[0] = prev_val[0];
        x[1] = prev_val[1];
        if (share_nu_alpha) {
            x[2] = prev_val[1];
        } else {
            x[2] = prev_val[2];
        }
        return;
    }

    x[0] = exp(par_out[0]);
    x[1] = exp(par_out[1]);
    if (share_nu_alpha) {
        x[2] = exp(par_out[1]);
    } else {
        x[2] = exp(par_out[2]);
    }

    if ((tau_sum < TAU_CHECK_VAL) && (!orig_share)) {
        x[2] = sqrt(DBL_MIN);
    }
    return;
}

void maximizePrecisionUnpaired(void * inputs, double * prev_val, double * x,
        bool share_nu_alpha) {
    // starting location
    precPtrsUnpaired * ex = static_cast<precPtrsUnpaired * >(inputs);
    int p = ex->p;
    int n_pep = ex->n_peptide;
    int n_trt = ex->n_trt;
    double TAU_CHECK_VAL = 1.0e-15;
    double tau_sum(0.0);
    bool orig_share = share_nu_alpha;

    for (int i = 0; i < n_trt; i++) {
        tau_sum += ex->tau[n_pep * i + p];
    }
    if ((tau_sum < TAU_CHECK_VAL) && !share_nu_alpha) {
        ex->share_nu_alpha = true;
        share_nu_alpha = true;
    }

    double par[3];
    double par_out[3];
    int n;

    par[0] = log(prev_val[0]);
    par[1] = log(prev_val[1]);
    if (share_nu_alpha) {
        par[2] = log(prev_val[1]);
        n = 2;
    } else {
        par[2] = log(prev_val[2]);
        n = 3;
    }

    // misc tuning parameters and reporting parameters
    int maxit(100), fncount, fail, trace(0);
    double Fmin, abstol(- INFINITY), intol(1.0e-4);
    double alpha(1.0), beta(.5), gamma(2.0);

    nmmin(n, par, par_out, &Fmin, precisionObjectiveFnUnpaired,
            &fail, abstol, intol, inputs, alpha, beta, gamma, trace,
            &fncount, maxit);

    if (fail > 0) {
        x[0] = prev_val[0];
        x[1] = prev_val[1];
        if (share_nu_alpha) {
            x[2] = prev_val[1];
        } else {
            x[2] = prev_val[2];
        }
        return;
    }

    x[0] = exp(par_out[0]);
    x[1] = exp(par_out[1]);
    if (share_nu_alpha) {
        x[2] = exp(par_out[1]);
    } else {
        x[2] = exp(par_out[2]);
    }

    if ((tau_sum < TAU_CHECK_VAL) && (!orig_share)) {
        x[2] = sqrt(DBL_MIN);
    }
    return;
}

double betaHyperObjectiveFn(int n, double * par, void * ex) {
    // m is location parameter, tau is log of precision parameter
    const double m(par[0]), tau(par[1]);
    // extract objective parameters
    double * input = static_cast<double *>(ex);
    const double beta1_sum(input[0]);
    const double beta1_sqr_sum(input[1]);
    const double P(input[2]);

    const double l1(input[3]);
    const double s1(input[4]);
    const double m_bar(input[5]);
    const double nu_m(input[6]);

    double out;
    out = - P * pnorm(m * exp(.5 * tau), 0, 1, /* lower_tail = */ 1,
            /* give_log = */ 1);
    out += - exp(tau) * (beta1_sqr_sum - 2.0 * m * beta1_sum + P * m * m) / 2.0;
    out += P * tau / 2.0;
    out += - (m - m_bar) * (m - m_bar) * nu_m / 2.0;
    out += -exp(tau) * l1 + (s1 - 1.0) * tau;
    return(- out);
}

void betaHyperObjectiveGr(int n, double * par, double * gr, void * ex) {
    // m is location parameter, tau is log of precision parameter
    const double m(par[0]), tau(par[1]);
    // extract objective parameters
    double * input = static_cast<double *>(ex);
    const double beta1_sum(input[0]);
    const double beta1_sqr_sum(input[1]);
    const double P(input[2]);

    const double l1(input[3]);
    const double s1(input[4]);
    const double m_bar(input[5]);
    const double nu_m(input[6]);

    double g_m, g_tau;
    const double log_inv_mills = dnorm(m * exp(.5 * tau), 0, 1, /* give_log */ 1) -
            pnorm(m * exp(.5 * tau), 0, 1, /* lower_tail */ 1, /* give_log */ 1);
    g_m = - P * exp(.5 * tau + log_inv_mills);
    g_m += exp(tau) * (beta1_sum - m * P);
    g_m += nu_m * (m_bar - m);

    g_tau = - P * m * .5 * exp(.5 * tau + log_inv_mills);
    g_tau += - exp(tau) * .5 * (beta1_sqr_sum - 2.0 * m * beta1_sum + P * m * m);
    g_tau += - exp(tau) * l1 + (s1 - 1.0) + P * .5;

    gr[0] = - g_m;
    gr[1] = - g_tau;
}

double gammaHyperObjectiveFn(int n, double * par, void * ex) {
    const double s(par[0]);
    if (s <= 0.0)
        return INFINITY;
    double * input = static_cast<double *>(ex);
    const double sum_log_x(input[0]);
    const double sum_x(input[1]);
    const double P(input[2]);
    const double l_s(input[3]);
    const double l_l(input[4]);

    double out = s * (-l_s + sum_log_x) - P * lgammafn(s) - P * s +
            P * s * log(P * s / (l_l + sum_x));
    return(- out);
}

void gammaHyperObjectiveGr(int n, double * par, double * gr, void * ex) {
    const double s(par[0]);
    double * input = static_cast<double *>(ex);
    const double sum_log_x(input[0]);
    const double sum_x(input[1]);
    const double P(input[2]);
    const double l_s(input[3]);
    const double l_l(input[4]);

    double g_s = -l_s + sum_log_x - P * digamma(s) +
            P * log(P * s / (l_l + sum_x));
    gr[0] = - g_s;
}

double precisionObjectiveFn(int n, double * par, void * ex) {
    // extract summary statistics and dimensions
    precPtrs * inputs = static_cast<precPtrs * >(ex);
    bool share = inputs->share_nu_alpha;
    const int n_subject = inputs->n_subject;
    const int n_isamp = inputs->n_isamp;
    const int n_isamp_max = inputs->n_isamp_max;
    const int n_peptide = inputs->n_peptide;
    const int p = inputs->p;

    const double n_rep = inputs->n_rep;
    const double beta0 = inputs->beta0;
    const double beta1 = inputs->beta1;
    const double s_a = inputs->s_a;
    const double s_e = inputs->s_e;
    const double l_a = inputs->l_a;
    const double l_e = inputs->l_e;

    const double l_nu_e(par[0]);
    const double l_nu_a0(par[1]);
    double l_nu_a1;
    if (share) {
        l_nu_a1 = par[1];
    } else {
        l_nu_a1 = par[2];
    }
    const double nu_e(exp(l_nu_e));
    const double nu_a0(exp(l_nu_a0));
    const double nu_a1(exp(l_nu_a1));


    double lik(0.0);
    for (int i = 0; i < n_subject; i++) {
        int k0, k1;
        double t0_w(0.0), t1_w(0.0), cross_w(0.0);
        double det_w(0.0), ss0_w(0.0), ss1_w(0.0);
        double det0_w, det1_w;
        double w0_g0, w1_g0, u0, q, bot;
        double w0_g1, w1_g1, u1, q0, q1, bot0, bot1;
        double tau1m, tau, ym0, ym1, yss0, yss1, mu0, mu1;

        k0 = n_peptide * i + p;
        tau = inputs->tau[k0];
        tau1m = 1.0 - tau;

        k0 = n_peptide * (2 * i) + p;
        k1 = n_peptide * (2 * i + 1) + p;

        ym0 = inputs->y_mean[k0];
        ym1 = inputs->y_mean[k1];
        yss0 = inputs->y_ss[k0];
        yss1 = inputs->y_ss[k1];
        mu0 = inputs->mu[2 * i];
        mu1 = inputs->mu[2 * i + 1];

        for (int r = 0; r < n_isamp; r++) {
            k0 = n_isamp_max * i + r;
            w0_g0 = inputs->w0_g0_samp[k0];
            w1_g0 = inputs->w1_g0_samp[k0];
            u0 = inputs->u0_g0_samp[k0];
            q = inputs->q0_g0_weights[k0];
            bot = n_rep * nu_e * (w0_g0 + w1_g0) + u0 * nu_a0;

            t0_w += w0_g0 * u0 * q / bot;
            t1_w += w1_g0 * u0 * q / bot;
            cross_w += w1_g0 * w0_g0 * q / bot;
            det_w += log(bot) * q;
            ss0_w += w0_g0 * q;
            ss1_w += w1_g0 * q;
        }

        lik += -.5 * det_w * tau1m;
        lik += -.5 * n_rep * nu_a0 * nu_e * t0_w * tau1m *
                pow2(ym0 - beta0 - mu0);
        lik += -.5 * n_rep * nu_a0 * nu_e * t1_w * tau1m *
                pow2(ym1 - beta0 - mu1);
        lik += -.5 * pow2(n_rep * nu_e) * pow2(ym0 - mu0 - ym1 + mu1) *
                tau1m * cross_w;
        lik += -.5 * yss0 * ss0_w * tau1m * nu_e;
        lik += -.5 * yss1 * ss1_w * tau1m * nu_e;

        t0_w = 0.0;
        t1_w = 0.0;
        det0_w = 0.0;
        det1_w = 0.0;
        ss0_w = 0.0;
        ss1_w = 0.0;

        for (int r = 0; r < n_isamp; r++) {
            k0 = n_isamp_max * i + r;
            w0_g1 = inputs->w0_g1_samp[k0];
            w1_g1 = inputs->w1_g1_samp[k0];
            u0 = inputs->u0_g1_samp[k0];
            u1 = inputs->u1_g1_samp[k0];
            q0 = inputs->q0_g1_weights[k0];
            q1 = inputs->q1_g1_weights[k0];
            bot0 = nu_e * n_rep * w0_g1 + u0 * nu_a0;
            bot1 = nu_e * n_rep * w1_g1 + u1 * nu_a1;

            t0_w += u0 * w0_g1 * q0 / bot0;
            t1_w += u1 * w1_g1 * q1 / bot1;
            det0_w += log(bot0) * q0;
            det1_w += log(bot1) * q1;
            ss0_w += w0_g1 * q0;
            ss1_w += w1_g1 * q1;
        }

        lik += -.5 * det0_w * tau;
        lik += -.5 * det1_w * tau;
        lik += -.5 * n_rep * nu_a0 * nu_e * t0_w * tau *
                pow2(ym0 - beta0 - mu0);
        lik += -.5 * n_rep * nu_a1 * nu_e * t1_w * tau *
                pow2(ym1 - beta0 - beta1 - mu1);
        lik += -.5 * yss0 * ss0_w * tau * nu_e;
        lik += -.5 * yss1 * ss1_w * tau * nu_e;
        lik += .5 * l_nu_a1 * tau;
    } // end parallel region
    lik += n_rep * n_subject * l_nu_e;
    lik += .5 * n_subject * l_nu_a0;

    //priors
    lik += (s_a - 1.0) * l_nu_a0 + (s_e - 1.0) * l_nu_e;
    lik += -l_a * nu_a0 - l_e * nu_e;
    if (!share){
        lik += -l_a * nu_a1 + (s_a - 1.0) * l_nu_a1;
    }
    // return negative value because optimizer minimizes;
    return -lik;
}

double precisionObjectiveFnUnpaired(int n, double * par, void * ex) {
    // extract summary statistics and dimensions
    precPtrsUnpaired * inputs = static_cast<precPtrsUnpaired * >(ex);
    bool share = inputs->share_nu_alpha;
    const int n_trt = inputs->n_trt;
    const int n_ctl = inputs->n_ctl;
    const int n_isamp = inputs->n_isamp;
    const int n_isamp_max = inputs->n_isamp_max;
    const int n_peptide = inputs->n_peptide;
    const int p = inputs->p;

    const double n_rep = inputs->n_rep;
    const double beta0 = inputs->beta0;
    const double beta1 = inputs->beta1;
    const double s_a = inputs->s_a;
    const double s_e = inputs->s_e;
    const double l_a = inputs->l_a;
    const double l_e = inputs->l_e;

    const double l_nu_e(par[0]);
    const double l_nu_a0(par[1]);
    double l_nu_a1;
    if (share) {
        l_nu_a1 = par[1];
    } else {
        l_nu_a1 = par[2];
    }
    const double nu_e(exp(l_nu_e));
    const double nu_a0(exp(l_nu_a0));
    const double nu_a1(exp(l_nu_a1));

    double lik(0.0);
    int k;
    for (int i = 0; i < n_trt; i++) {
        double t0_w(0.0), t1_w(0.0);
        double ss0_w(0.0), ss1_w(0.0);
        double det0_w(0.0), det1_w(0.0);
        double w_g0, w_g1, u_g0, u_g1, q_g0, q_g1;
        double bot0, bot1;
        double tau1m, tau, ym, yss, mu;

        k = n_peptide * i + p;
        tau = inputs->tau[k];
        tau1m = 1.0 - tau;

        ym = inputs->y_mean_trt[k];
        yss = inputs->y_ss_trt[k];
        mu = inputs->mu[i + n_ctl];

        for (int r = 0; r < n_isamp; r++) {
            k = n_isamp_max * i + r;
            w_g0 = inputs->w_trt_g0[k];
            w_g1 = inputs->w_trt_g1[k];
            u_g0 = inputs->u_trt_g0[k];
            u_g1 = inputs->u_trt_g1[k];
            q_g0 = inputs->q_trt_g0[k];
            q_g1 = inputs->q_trt_g1[k];
            bot0 = nu_e * n_rep * w_g0 + u_g0 * nu_a0;
            bot1 = nu_e * n_rep * w_g1 + u_g1 * nu_a1;

            t0_w += u_g0 * w_g0 * q_g0 / bot0;
            t1_w += u_g1 * w_g1 * q_g1 / bot1;
            det0_w += log(bot0) * q_g0;
            det1_w += log(bot1) * q_g1;
            ss0_w += w_g0 * q_g0;
            ss1_w += w_g1 * q_g1;
        }

        lik += -.5 * det0_w * tau1m;
        lik += -.5 * det1_w * tau;
        lik += -.5 * n_rep * nu_a0 * nu_e * t0_w * tau1m *
                pow2(ym - beta0 - mu);
        lik += -.5 * n_rep * nu_a1 * nu_e * t1_w * tau *
                pow2(ym - beta0 - beta1 - mu);
        lik += -.5 * yss * ss0_w * tau1m * nu_e;
        lik += -.5 * yss * ss1_w * tau * nu_e;
        lik += .5 * l_nu_a1 * tau;
        lik += .5 * l_nu_a0 * tau1m;
    }
    lik += .5 * n_rep * n_trt * l_nu_e;

    double t, det, ss_w;
    double w, u, q;
    double bot;
    double ym, yss, mu;
    for (int i = 0; i < n_ctl; i++) {
        t = 0.0;
        det = 0.0;
        ss_w = 0.0;
        k = n_peptide * i + p;
        ym = inputs->y_mean_ctl[k];
        yss = inputs->y_ss_ctl[k];
        mu = inputs->mu[i];

        for (int r = 0; r < n_isamp; r++) {
            k = n_isamp_max * i + r;
            w = inputs->w_ctl[k];
            u = inputs->u_ctl[k];
            q = inputs->q_ctl[k];
            bot = nu_e * n_rep * w + u * nu_a0;

            t += u * w * q / bot;
            det += log(bot) * q;
            ss_w += w * q;
        }
        lik += -.5 * det;
        lik += -.5 * n_rep * nu_a0 * nu_e * t *
                pow2(ym - beta0 - mu);
        lik += -.5 * yss * ss_w * nu_e;
    }
    lik += .5 * n_rep * n_ctl * l_nu_e;
    lik += .5 * n_ctl * l_nu_a0;

    //priors
    lik += (s_a - 1.0) * l_nu_a0 + (s_e - 1.0) * l_nu_e;
    lik += -l_a * nu_a0 - l_e * nu_e;
    if (!share){
        lik += -l_a * nu_a1 + (s_a - 1.0) * l_nu_a1;
    }
    // return negative value because optimizer minimizes;
    return -lik;
}
