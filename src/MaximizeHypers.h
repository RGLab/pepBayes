/*
 * MaximizeHypers.h
 *
 *  Created on: Jun 27, 2014
 *      Author: gimholte
 */

#ifndef MAXIMIZEHYPERS_H_
#define MAXIMIZEHYPERS_H_

struct precPtrsUnpaired {
    double * y_mean_trt;
    double * y_mean_ctl;
    double * y_ss_trt;
    double * y_ss_ctl;
    double * mu;
    double * tau;

    double * q_trt_g0;
    double * q_trt_g1;
    double * q_ctl;

    double * w_trt_g0;
    double * w_trt_g1;
    double * w_ctl;

    double * u_trt_g0;
    double * u_trt_g1;
    double * u_ctl;

    double s_e;
    double l_e;
    double s_a;
    double l_a;
    double beta0;
    double beta1;

    double n_rep;
    int n_isamp;
    int n_isamp_max;
    int p;
    int n_trt;
    int n_ctl;
    int n_peptide;
    int n_threads;

    bool share_nu_alpha;
};

struct precPtrs {
    double * y_mean;
    double * y_ss;
    double * mu;
    double * tau;

    double * q0_g0_weights;
    double * q0_g1_weights;
    double * q1_g1_weights;

    double * w0_g0_samp;
    double * w0_g1_samp;
    double * w1_g0_samp;
    double * w1_g1_samp;

    double * u0_g0_samp;
    double * u0_g1_samp;
    double * u1_g1_samp;

    double s_e;
    double l_e;
    double s_a;
    double l_a;
    double beta0;
    double beta1;

    double n_rep;
    int n_isamp;
    int n_isamp_max;
    int p;
    int n_subject;
    int n_peptide;
    int n_threads;

    bool share_nu_alpha;
};

void maximizePrecisionUnpaired(void * inputs, double * prev_val, double * x,
        bool share_nu_alpha);
void maximizePrecision(void * inputs, double * prev_val, double * x,
        bool share_nu_alpha);
void maximizeGammaHypers(const double * data, int data_len, const double l_s,
        const double l_lambda, double * x);
void maximizeTnormHypers(const double * data, int data_len, const double m_mean,
        const double m_prec, const double nu_rate, const double nu_shape,
        double * x);

double betaHyperObjectiveFn(int n, double * par, void * ex);
void betaHyperObjectiveGr(int n, double * par, double * gr, void * ex);

double gammaHyperObjectiveFn(int n, double * par, void * ex);
void gammaHyperObjectiveGr(int n, double * par, double * gr, void * ex);

double precisionObjectiveFn(int n, double * par, void * ex);
double precisionObjectiveFnUnpaired(int n, double * par, void * ex);

#endif /* MAXIMIZEHYPERS_H_ */
