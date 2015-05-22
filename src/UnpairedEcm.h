/*
 * UnpairedEcm.h
 *
 *  Created on: Jun 12, 2014
 *      Author: gimholte
 */

#ifndef UNPAIREDECM_H_
#define UNPAIREDECM_H_

class UnpairedEcm: public EcmAlgorithm {
public:
    UnpairedEcm(const Rcpp::List & par_list, const  Rcpp::List & data_list,
            const Rcpp::List & iter_list);
    UnpairedEcm();
    virtual ~UnpairedEcm() {};
private:
    // model data
    arma::mat y_mean_trt, y_mean_ctl;
    arma::mat y_ss_trt, y_ss_ctl;

    // sample sizes
    int n_position, n_peptide, n_trt, n_ctl, n_slide;
    arma::vec n_rep;
    arma::ivec pos_start, pos_count, pos;

    // location parameters;
    arma::vec beta0, beta1, mu;

    // mixture parameters
    arma::mat tau;
    arma::vec a, b, omega;

    // precision parameters;
    arma::vec nu_alpha0, nu_alpha1, nu_eps;

    // estimated precision hyperparameters;
    double nu_beta0, nu_beta1;
    double s_alpha, s_eps;
    double lambda_alpha, lambda_eps;

    // estimated location hyperparameters;
    double m_beta0, m_beta1;

    // fixed location hyperparameters;
    const double m_beta0_prior_mean, m_beta0_prior_prec;
    const double m_beta1_prior_mean, m_beta1_prior_prec;
    const double nu_beta0_shape, nu_beta0_rate;
    const double nu_beta1_shape, nu_beta1_rate;
    const double s_alpha_rate, s_eps_rate;
    const double lambda_alpha_rate, lambda_eps_rate;
    const double lambda_a, lambda_b;

    // t-distribution parameters
    arma::cube w_ctl, w_trt_g0, w_trt_g1;
    arma::cube u_ctl, u_trt_g0, u_trt_g1;
    arma::mat shrink_wu_ctl_avg, shrink_wu_g0_avg, shrink_wu_g1_avg;
    arma::mat mean_shrink_g0, mean_shrink_g1, mean_shrink_ctl;
    arma::mat aprec_g0, aprec_g1, aprec_ctl;

    // marginal density and importance sampling parameters
    int n_isamp, n_rule, n_isamp_max;
    arma::ivec n_isamp_schedule;
    arma::cube q_ctl, q_trt_g0, q_trt_g1;
    arma::vec gauss_rule_x, gauss_rule_w;
    std::vector<arma::vec> gauss_work_x, gauss_work_gx;
    std::vector<arma::vec> alpha_is_work, density_inst_work, density_joint_work;
    std::vector<RngStream> rng;
    double nu;

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

    bool update_omega, update_beta_hypers;
    bool update_mu, update_beta0, update_beta1;
    bool update_precision, update_precision_hypers;
    bool update_weights;
    bool share_nu_alpha;

    // methods
    void iterate(const int & iter_idx);
    void collectIteration(const int & iter_idx);
    const Rcpp::List gatherOutput();
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
    void prepMuTrt(const int & i);
    void prepMuCtl(const int & i);
    void updateWeightSummaries(const int & p);
    void setImportanceSamplingSchedule();

    void sampleWeightsTrt(const int & p, const int & i, const int & th_id,
            const double & delta_g0, const double & delta_g1,
            const double & offset, arma::vec & alpha_is_work,
            arma::vec & density_inst_work,
            arma::vec & density_joint_work);

    void sampleWeightsCtl(const int & p, const int & i, const int & th_id,
            const double & delta, const double & offset,
            arma::vec & alpha_is_work,
            arma::vec & density_inst_work,
            arma::vec & density_joint_work);
};

#endif /* UNPAIREDECM_H_ */
