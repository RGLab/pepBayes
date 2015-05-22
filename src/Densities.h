/*
 * Densities.h
 *
 *  Created on: Jun 12, 2014
 *      Author: gimholte
 */

#ifndef DENSITIES_H_
#define DENSITIES_H_

#include <math.h>

template <typename derived>
class LogDensity {
public:
    double pdf(const double & x) {
        return static_cast<derived*>(this)->pdf(x);
    }

    double pdfDeriv(const double & x) {
        return static_cast<derived*>(this)->pdfDeriv(x);
    }

    void setParameters(double const * const args){
        static_cast<derived*>(this)->setParameters(args);
    }

    void printParameters() {
        static_cast<derived*>(this)->printParameters();
    }
};

class ShapeConditional : public LogDensity<ShapeConditional> {
    double p, log_of_sum, sum_of_log_nu, prior_s_rate;
public:
    void setParameters(double const * const args) {
        p = args[0];
        log_of_sum = args[1];
        sum_of_log_nu = args[2];
        prior_s_rate = args[3];
    }
    double pdf(const double & x) const {
        return R::lgammafn(p * x + 1.0) -
                p * R::lgammafn(x) - x * (prior_s_rate - sum_of_log_nu) -
                p * x * log_of_sum;
    }

    void printParameters() const {
        Rcpp::Rcout << " " << p << " " << log_of_sum << " " <<
                sum_of_log_nu << std::endl;
    }

    double pdfDeriv(const double & x) const {
        return p * R::digamma(p * x + 1.0) - p * R::digamma(x) -
                (prior_s_rate - sum_of_log_nu) - p * log_of_sum;
    }
};

class BetaHyperConditional : public LogDensity<BetaHyperConditional> {
    double q, twin_par, log_sum;
public:
    void setParameters(double const * const args) {
        q = args[0];
        twin_par = args[1];
        log_sum = args[2];
    }
    double pdf(const double & x) const {
        return - q * R::lbeta(x, twin_par) - x * log_sum;
    }
    void printParameters() const {
        Rcpp::Rcout << " " << q << " " << log_sum  << " " << twin_par << std::endl;
    }
    double pdfDeriv(const double & x) const {
        return - log_sum + q * R::digamma(x + twin_par) - q * R::digamma(x);
    }
};

#endif /* DENSITIES_H_ */
