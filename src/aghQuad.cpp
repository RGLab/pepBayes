/*
 * aghQuad.cpp
 *
 *  Created on: Nov 11, 2014
 *      Author: gimholte
 */

#include "RcppArmadillo.h"
#include "aghQuad.h"

RcppExport SEXP testAgh(SEXP r_parameter_list, SEXP r_data_list,
        SEXP r_gauss_rule_x, SEXP r_gauss_rule_w) {
    Rcpp::List par_list = Rcpp::as<Rcpp::List>(r_parameter_list);
    Rcpp::List data_list = Rcpp::as<Rcpp::List>(r_data_list);

    const double nu_eps = par_list["nu_eps"];
    const double nu_alpha0 = par_list["nu_alpha0"];
    const double nu_alpha1 = par_list["nu_alpha1"];
    const double beta0 = par_list["beta0"];
    const double mu0 = par_list["mu0"];
    const double beta1 = par_list["beta1"];
    const double mu1 = par_list["mu1"];
    const double nu = par_list["nu"];

    const double y_ss0 = data_list["y_ss0"];
    const double y_ss1 = data_list["y_ss1"];
    const double y_mean0 = data_list["y_mean0"];
    const double y_mean1 = data_list["y_mean1"];
    const double n_rep = data_list["n_rep"];

    const double delta_0 = y_mean0 - beta0 - mu0;
    const double delta_1_0 = y_mean1 - beta0 - mu1;
    const double delta_1_1 = y_mean1 - beta0 - beta1 - mu1;

    const arma::vec gauss_rule_x = Rcpp::as<arma::vec>(r_gauss_rule_x);
    const arma::vec gauss_rule_w = Rcpp::as<arma::vec>(r_gauss_rule_w);
    const int n_rule = (int) gauss_rule_w.n_elem;
    arma::vec gauss_work_x(n_rule);
    arma::vec gauss_work_gx(n_rule);

    arma::vec out_density(2);
    const double prec0(n_rep * nu_eps + nu_alpha0);
    const double prec1(n_rep * nu_eps + nu_alpha1);
    const double shrink0(n_rep * nu_eps / prec0);
    const double shrink1(n_rep * nu_eps / prec1);
    double out[3];
    out_density(0) = densityGamma0(/* offset0 = */ y_ss0,
            /* delta0 = */ delta_0, /* offset1 = */ y_ss1,
            /* delta1 = */ delta_1_0, nu, nu_eps, nu_alpha0,
            n_rep, gauss_rule_x, gauss_rule_w,
            gauss_work_x, gauss_work_gx);
    arma::vec gauss_work_gx0(gauss_work_gx);

    densityGamma1(out, /* offset0 = */ y_ss0,
            /* delta0 = */ delta_0, /* offset1 = */ y_ss1,
            /* delta1 = */ delta_1_1, nu, nu_eps, nu_alpha0,
            nu_alpha1, n_rep, delta_0 * shrink0, delta_1_1 * shrink1,
            1.0 / sqrt(prec0), 1.0 / sqrt(prec1), gauss_rule_x, gauss_rule_w,
            gauss_work_x, gauss_work_gx);
    out_density(1) = out[2];
    return Rcpp::List::create(Rcpp::Named("dens_value") = out_density,
            Rcpp::Named("work0") = gauss_work_gx0,
            Rcpp::Named("work1") = gauss_work_gx);
}

double aghQuad(const arma::vec & gauss_rule_w,
        const arma::vec & gauss_rule_x,
        const arma::vec & gauss_work_gx,
        const double & sigma) {
    double I = 0.0;
    const int n = (int) gauss_rule_w.n_elem;
    for (int i = 0; i < n; i++) {
        I += exp(gauss_rule_x(i) *
                gauss_rule_x(i) + gauss_work_gx(i) + gauss_rule_w(i));
    }
    return(M_SQRT2 * sigma * I);
}

// evaluate the marginal density of an observation set
// given gamma = 1 (i.e, random effects and weights integrated out)
void densityGamma1(double * out, const double & offset0, const double & delta0,
        const double & offset1, const double & delta1, const int & nu,
        const double & nu_eps, const double & nu_alpha0, const double & nu_alpha1,
        const double & n_rep,
        const double & m0, const double & m1,
        const double & sigma0, const double & sigma1,
        const arma::vec & gauss_rule_x, const arma::vec & gauss_rule_w,
        arma::vec & gauss_work_x, arma::vec & gauss_work_gx) {
    // compute the joint density of the data and gamma_ip = 0
    // for all subjects within a single peptide.
    // Uses Gauss-Hermite quadrature

    // translate eval points and then evaluate joint at the
    // translated points
    densityAlphaGamma1 f1(nu, n_rep);
    gauss_work_x = m0 + M_SQRT2 * gauss_rule_x * sigma0;
    f1.calc(gauss_work_gx, gauss_work_x,
            offset0, delta0, nu, nu_eps, nu_alpha0, n_rep);
    out[0] = aghQuad(gauss_rule_w, gauss_rule_x, gauss_work_gx, sigma0);

    gauss_work_x = m1 + M_SQRT2 * gauss_rule_x * sigma1;
    f1.calc(gauss_work_gx, gauss_work_x,
            offset1, delta1, nu, nu_eps, nu_alpha1, n_rep);
    out[1] = aghQuad(gauss_rule_w, gauss_rule_x, gauss_work_gx, sigma1);
    out[2] = out[0] * out[1];
}


// evaluate the marginal density of an observation set
// given gamma = 0 (i.e, random effects and weights integrated out)
double densityGamma0(const double & offset0, const double & delta0,
        const double & offset1, const double & delta1, const int & nu,
        const double & nu_eps, const double & nu_alpha,
        const double & n_rep,
        const arma::vec & gauss_rule_x, const arma::vec & gauss_rule_w,
        arma::vec & gauss_work_x, arma::vec & gauss_work_gx) {
    // compute the joint density of the data and gamma_ip = 0
    // for all subjects within a single peptide.
    // Uses Gauss-Hermite quadrature
    const double nu_star = 2.0 * n_rep * nu_eps + nu_alpha;
    const double sigma = sqrt(nu / ((nu - 1) * nu_star));
    const double m_star = n_rep * nu_eps * (delta0 + delta1) / nu_star;

    // translate eval points and then evaluate joint at the
    // translated points
    gauss_work_x = m_star + M_SQRT2 * gauss_rule_x * sigma;
    densityAlphaGamma0 f0(nu, n_rep);
    f0.calc(gauss_work_gx, gauss_work_x,
            offset0, delta0, offset1, delta1, nu, nu_eps, nu_alpha,
            n_rep);

    return aghQuad(gauss_rule_w,
            gauss_rule_x, gauss_work_gx, sigma);
}

void densityAlphaGamma1::calc(arma::vec & out, const arma::vec & alpha,
        const double & offset, const double & delta,
        const int & nu, const double & nu_eps, const double & nu_alpha,
        const double & n_rep) {
    const double scale_data =
            (n_rep / 2.0) * log(nu_eps / (M_PI * nu));
    const double scale_ranef = .5 * log(nu_alpha / (M_PI * nu));
    double den;
    for (unsigned j = 0; j < alpha.n_elem; j++) {
        den = -.5 * (n_rep + nu) * log1p(nu_eps *
                (offset + n_rep * pow(delta - alpha[j], 2)) / nu);
        den += -.5 * (1.0 + nu) *
                log1p(nu_alpha * alpha[j] * alpha[j] / nu);
        den += log_gamma_ratio_data + log_gamma_ratio_ranef;
        den += scale_data + scale_ranef;
        out[j] = den;
    }
    return;
}

void densityAlphaGamma0::calc(arma::vec & out, const arma::vec & alpha,
        const double & offset0, const double & delta0,
        const double & offset1, const double & delta1,
        const int & nu, const double & nu_eps, const double & nu_alpha,
        const double & n_rep) {
    const double scale_data = (n_rep / 2.0) *
            log(nu_eps / (M_PI * nu));
    const double scale_ranef = .5 * log(nu_alpha / (M_PI * nu));
    double den;
    for (unsigned j = 0; j < alpha.n_elem; j++) {
        den = -.5 * (n_rep + nu) * log1p(nu_eps *
                (offset0 + n_rep * pow(delta0 - alpha[j], 2)) / nu);
        den += -.5 * (n_rep + nu) * log1p(nu_eps *
                (offset1 + n_rep * pow(delta1 - alpha[j], 2)) / nu);
        den += -.5 * (1.0 + nu) * log1p(nu_alpha * alpha[j] * alpha[j] / nu);
        den += 2.0 * log_gamma_ratio_data + log_gamma_ratio_ranef;
        den += 2.0 * scale_data + scale_ranef;
        out[j] = den;
    }
    return;
}
