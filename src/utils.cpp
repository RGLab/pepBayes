/*
 * utils.cpp
 *
 *  Created on: Sep 27, 2013
 *      Author: hoblitz
 */
#include "RngStream.h"
#include <Rcpp.h>
#include "utils.h"

double logspaceAdd(const double loga, const double logb) {
    if (!R_FINITE(loga))
        return logb;
    if (loga > logb)
        return logspaceAdd(logb, loga);
    return logb + log1p(exp(loga - logb));
}

double logspaceSubtract(const double loga, const double logb) {
    if (logb >= loga) {
        Rcpp::stop("Must have b < a for logspace subtraction.");
    }
    if (!R_FINITE(logb) && logb < 0.0)
        return loga;
    return loga + log1p(-exp(logb - loga));
}

bool metropolisHastings(const double & cur_value, const double & cur_lik,
        const double & cur_to_prop_trans, const double & prop_value,
        const double & prop_lik, const double & prop_to_cur_trans,
        RngStream rng) {
    const double loga = prop_lik + prop_to_cur_trans - cur_lik - cur_to_prop_trans;
    if (loga >= log(RngStream_RandU01(rng)))
        return true;
    else
        return false;
}

double logFromLogit(const double logit_x) {
    return logit_x - R::log1pexp(logit_x);
}

double log1mFromLogit(const double logit_x) {
    return - R::log1pexp(logit_x);
}

double pow2(const double & x) {
    return x * x;
}
