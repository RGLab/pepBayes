/*
 * PairedEcm.h
 *
 *  Created on: Nov 3, 2014
 *      Author: gimholte
 */

#ifndef PAIREDECM_H_
#define PAIREDECM_H_

class PairedEcm: public EcmAlgorithm {
public:
    PairedEcm(const Rcpp::List & par_list, const Rcpp::List & data_list,
            const Rcpp::List & iter_list);
    virtual ~PairedEcm() {};
private:
    // model data
    arma::mat y_mean, y_ss;
    int n_position, n_slide, n_subject, n_peptide;
    arma::vec n_rep;
    arma::ivec pos_start, pos_count, pos;

    // location parameters;
    arma::vec beta0, beta1, mu;

    // estimated location hyperparameters;
    double m_beta0, m_beta1;

    // fixed hyperparameters;
    const double m_beta0_prior_mean, m_beta0_prior_prec;
    const double m_beta1_prior_mean, m_beta1_prior_prec;
    const double s_alpha_rate, s_eps_rate;
    const double lambda_alpha_rate, lambda_eps_rate;
    const double nu_beta0_shape, nu_beta0_rate;
    const double nu_beta1_shape, nu_beta1_rate;
    const double lambda_a, lambda_b;

    // iteration controls for fixing parameters
    bool update_omega;
    bool update_mu, update_beta0, update_beta1, update_beta_hypers;
    bool update_precision, update_precision_hypers;
    bool update_weights;
    bool share_nu_alpha;

    // precision parameters;
    arma::vec nu_alpha0, nu_alpha1, nu_eps;

    // estimated precision hyperparameters;
    double nu_beta0, nu_beta1;
    double s_alpha, s_eps;
    double lambda_alpha, lambda_eps;

    // t-distribution parameters
    arma::cube w0_g0_samp, w0_g1_samp, w1_g0_samp, w1_g1_samp;
    arma::cube u0_g0_samp, u0_g1_samp, u1_g1_samp;
    arma::mat shrink_w0_g0_avg, shrink_w1_g0_avg;
    arma::mat shrink_w0_g1_avg, shrink_w1_g1_avg;
    arma::mat shrink_cross_g0_avg;
    const double nu;

    // marginal density and importance sampling parameters
    int n_isamp, n_rule, n_isamp_max;
    arma::ivec n_isamp_schedule;
    arma::cube q0_g0_weights, q0_g1_weights, q1_g1_weights;
    arma::vec gauss_rule_x, gauss_rule_w;
    std::vector<arma::vec> gauss_work_x, gauss_work_gx;
    std::vector<arma::vec> alpha_is_work, density_inst_work, density_joint_work;
    std::vector<RngStream> rng;

    // mixture parameters
    arma::mat tau;
    arma::vec a, b, omega;

    // mu update helpers
    arma::vec mu_mean, mu_temp, mu_prec_inner_diag;
    arma::mat Q, mu_omega;
    arma::mat mu_prec_inner;

    // likelihood and trace storage;
    arma::mat hypers_tr;
    arma::mat beta0_tr;
    arma::mat beta1_tr;
    arma::mat mu_tr;
    arma::mat nu_eps_tr;
    arma::mat nu_alpha0_tr;
    arma::mat nu_alpha1_tr;
    arma::mat omega_tr;

    // methods
    void iterate(const int & iter_idx);
    void collectIteration(const int & iter_idx);
    const Rcpp::List gatherOutput();

    void sampleWeightsGamma0(const int & p,
            const int & i, const int & th_id,
            const double & delta0, const double & delta1_g0,
            const double & offset0, const double & offset1,
            arma::vec & alpha_is_work, arma::vec & density_inst_work,
            arma::vec & density_joint_work);

    void sampleWeightsGamma1(const int & p,
            const int & i, const int & th_id,
            const double & delta0, const double & delta1_g0,
            const double & offset0, const double & offset1,
            arma::vec & alpha_is_work, arma::vec & density_inst_work,
                        arma::vec & density_joint_work);

    bool checkInput();
    void updateBeta0(const int & p);
    void updateBeta1(const int & p);
    void updateM0();
    void updateNuBeta0();
    void updatePrecision(const int & p);
    void updateBeta1Hypers();
    void updatePrecHypers();
    void updateOmega(const int & p);
    void eStep(const int & p, const int & th_id);
    void updateMu();
    void prepMu(const int & i);
    void updateWeightSummaries(const int & p);
    void setImportanceSamplingSchedule();
};

#endif /* TESTPAIREDECM_H_ */
