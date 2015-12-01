/*
 * UnpairedMcmc.h
 *
 *  Created on: Nov 21, 2014
 *      Author: gimholte
 */

#ifndef UNPAIREDMCMC_H_
#define UNPAIREDMCMC_H_

class UnpairedMcmc: public MarkovChain {
    arma::mat y_mean_trt;
    arma::mat y_mean_control;
    arma::mat y_ss_trt;
    arma::mat y_ss_control;
    arma::ivec n_rep;

    arma::ivec pos;
    arma::ivec pos_start;
    arma::ivec pos_count;

    int n_position;
    int n_peptide;
    int n_slide;
    int n_control;
    int n_trt;

    arma::mat gamma;
    arma::mat gamma_prob;
    arma::vec omega;
    arma::vec a;
    arma::vec b;

    arma::vec beta0;
    arma::vec beta1;
    arma::vec mu;
    arma::mat alpha_trt;
    arma::mat alpha_ctl;

    arma::vec nu_eps;
    arma::vec nu_alpha0;
    arma::vec nu_alpha1;

    double nu_beta0;
    double m_beta0;
    double nu_beta1;
    double m_beta1;
    double s_alpha;
    double lambda_alpha;
    double s_eps;
    double lambda_eps;

    const double nu_beta0_shape;
    const double nu_beta0_rate;
    const double m_beta0_prior_mean;
    const double m_beta0_prior_prec;
    const double nu_beta1_shape;
    const double nu_beta1_rate;
    const double m_beta1_prior_mean;
    const double m_beta1_prior_prec;
    const double lambda_alpha_rate;
    const double lambda_eps_rate;
    const double s_alpha_rate;
    const double s_eps_rate;
    const double lambda_a;
    const double lambda_b;

    bool update_omega, update_omega_hypers;
    bool update_mu, update_beta0, update_beta1, update_beta_hypers;
    bool update_precision, update_precision_hypers;
    bool update_alpha, update_weights, update_gamma, update_dof;

    arma::mat w_trt;
    arma::mat w_ctl;
    arma::mat u_trt;
    arma::mat u_ctl;
    arma::mat mahala_dist_trt, mahala_dist_ctl;
    double nu_err, nu_err_transform;
    double nu_re, nu_re_transform;

    arma::vec mu_mean;
    arma::vec mu_star_mean;
    arma::vec mu_star;
    arma::vec mu_prec_inner_diag;
    arma::vec mu_temp;
    arma::mat Q, mu_omega;

    MHTuner m_beta1_tuner, nu_beta1_tuner, nu_re_tuner, nu_err_tuner;
    std::vector<AdaptiveGibbsSampler<BetaHyperConditional> > a_sampler;
    std::vector<AdaptiveGibbsSampler<BetaHyperConditional> > b_sampler;

    AdaptiveGibbsSampler<ShapeConditional> s_eps_sampler;
    AdaptiveGibbsSampler<ShapeConditional> s_alpha_sampler;
    std::vector<RngStream> rng;

    /* trace storage */
    arma::mat beta0_trace;
    arma::mat beta1_trace;
    arma::mat mu_trace;
    arma::mat a_trace;
    arma::mat b_trace;
    arma::mat hypers_trace;
    arma::mat nu_alpha0_trace;
    arma::mat nu_alpha1_trace;
    arma::mat nu_eps_trace;
    arma::mat omega_trace;
    arma::mat ppa;

    /* private methods */
    void computeMuMean_Trt(const int & i);
    void computeMuMean_Control(const int & i);
    void updateMu(RngStream & rng);
    void updateGammaAlphaWeights_Trt(const int & i, RngStream & rng);
    void updateAlphaWeights_Ctl(const int & i, RngStream & rng);
    void updateM0(RngStream & rng);
    void updateBeta0(const int & p, RngStream & rng);
    void updateBeta1(const int & p, RngStream & rng);
    void updateAlphaHypers(RngStream & rng);
    void updateEpsilonHypers(RngStream & rng);
    void updateNuBeta(RngStream & rng);
    void updateNuEps(const int & p, RngStream & rng);
    void updateNuAlpha(const int & p, RngStream & rng);
    void updateOmega(const int & p, RngStream & rng);
    void updateM1(RngStream & rng);
    void updateNuBeta1(RngStream & rng);
    void updateBetaHypers(const int & q, RngStream & rng);
    void updateNuErr(RngStream & rng);
    void updateNuRe(RngStream & rng);
    void computeMahalaDist(const int & i);
    void mahalaDistTrt(const int & i);
    void mahalaDistCtl(const int & i);
    double dfReDensity(const double x);
    double dfErrDensity(const double x);


public:
    UnpairedMcmc(const Rcpp::List & prior_pars, const Rcpp::List & data_list,
            const Rcpp::List & chain_pars);
    void iterate();
    void collectIteration(const int & sample_idx);
    const Rcpp::List chainOutput();
};


#endif /* TESTUNPAIREDMCMC_H_ */
