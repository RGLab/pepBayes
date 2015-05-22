/*
 * utils.h
 *
 *  Created on: Sep 27, 2013
 *      Author: hoblitz
 */

#ifndef UTILS_H_
#define UTILS_H_

double logspaceAdd(const double loga, const double logb);
double logspaceSubtract(const double loga, const double logb);
double metropolisHastings(const double & cur_value, const double & cur_lik,
        const double & cur_to_prop_trans, const double & prop_value,
        const double & prop_lik, const double & prop_to_cur_trans,
        RngStream rng);
double logFromLogit(const double logit_x);
double log1mFromLogit(const double logit_x);
double pow2(const double & x);
#endif /* UTILS_H_ */
