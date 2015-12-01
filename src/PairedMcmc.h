/*
 * PairedMcmc.h
 *
 *  Created on: Oct 23, 2014
 *      Author: gimholte
 */

#ifndef PAIREDMCMC_H_
#define PAIREDMCMC_H_

class PairedMcmc: public MarkovChain {
    arma::mat y_mean;
    arma::mat y_ss;
    arma::ivec n_rep;

    arma::ivec pos;
    arma::ivec pos_start;
    arma::ivec pos_count;

    int n_position;
    int n_peptide;
    int n_slide;
    int n_subject;

    arma::vec beta0;
    arma::vec beta1;
    arma::vec mu;
    arma::mat alpha0;
    arma::mat alpha1;

    arma::mat gamma;
    arma::mat gamma_prob;
    arma::vec omega;
    arma::vec a;
    arma::vec b;

    arma::vec nu_alpha0;
    arma::vec nu_alpha1;
    arma::vec nu_eps;

    double nu_beta0;
    double m_beta0;
    double nu_beta1;
    double m_beta1;
    double s_alpha;
    double lambda_alpha;
    double s_eps;
    double lambda_eps;

    // fixed hyperparameters
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
    bool update_weights, update_gamma, update_alpha, update_dof;

    bool share_nu_alpha;

    // t-distribution dof estimation helpers
    double nu_err, nu_err_transform;
    double nu_re, nu_re_transform;
    MHTuner nu_err_tuner;
    MHTuner nu_re_tuner;
    arma::mat w;
    arma::mat u;
    arma::mat mahala_dist;

    // helpers for sampling mu
    arma::vec mu_mean;
    arma::vec mu_star_mean;
    arma::vec mu_star;
    arma::vec mu_prec_inner_diag;
    arma::vec mu_temp;
    arma::mat Q, mu_omega, Qt;

    // helpers for sampling misc parameters
    MHTuner m_beta1_tuner, nu_beta1_tuner;
    std::vector<AdaptiveGibbsSampler<BetaHyperConditional> > a_sampler;
    std::vector<AdaptiveGibbsSampler<BetaHyperConditional> > b_sampler;
    AdaptiveGibbsSampler<ShapeConditional> s_eps_sampler;
    AdaptiveGibbsSampler<ShapeConditional> s_alpha_sampler;
    std::vector<RngStream> rng;

    /* trace storage */
    arma::mat beta0_trace;
    arma::mat beta1_trace;
    arma::mat mu_trace;

    arma::mat omega_trace;
    arma::mat a_trace;
    arma::mat b_trace;

    arma::mat nu_eps_trace;
    arma::mat nu_alpha0_trace;
    arma::mat nu_alpha1_trace;

    arma::mat hypers_trace;
    arma::mat ppa;
    arma::mat fitted_response;

    /* methods */
    void computeMuMean(const int & i);
    void updateMu(RngStream & rng);
    void updateM0(RngStream & rng);
    void updateNuBeta0(RngStream & rng);
    void updateBeta0(const int & p, RngStream & rng);
    void updateBeta1(const int & p, RngStream & rng);
    void updateAlphaHypers(RngStream & rng);
    void updateEpsilonHypers(RngStream & rng);
    void updateNuEps(const int & p, RngStream & rng);
    void updateNuAlpha(const int & p, RngStream & rng);
    void updateOmega(const int & p, RngStream & rng);
    void updateM1(RngStream & rng);
    void updateNuBeta1(RngStream & rng);
    void updateBetaHypers(const int & q, RngStream & rng);
    void updateGammaAlpha(const int & i, RngStream & rng);
    void computeMahalaDist(const int & i);
    void updateWeights(const int & i, RngStream & rng);
    void updateNuErr(RngStream & rng);
    void updateNuRe(RngStream & rng);
    double dfReDensity(const double x);
    double dfJointDensity(const double x);


public:
    PairedMcmc(const Rcpp::List & par_list, const Rcpp::List & data_list,
            const Rcpp::List & chain_list);
    void iterate();
    void collectIteration(const int & sample_idx);
    const Rcpp::List chainOutput();
};


#endif /* PAIREDMCMC_H_ */
