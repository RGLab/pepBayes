/*
 * norm_gamma_generation.h
 *
 *  Created on: Oct 6, 2011
 *      Author: hoblitz
 */

#ifndef NORM_GAMMA_GENERATION_H_
#define NORM_GAMMA_GENERATION_H_

#endif /* NORM_GAMMA_GENERATION_H_ */

double RngStream_N01 (const RngStream r);
double RngStream_GA1 (const double a, RngStream r);
double RngStream_Beta (const double a, const double b, RngStream r);
double RngStream_LogitBeta(const double a, const double b, RngStream rng);
double RngStream_TruncNorm(const double & mean, const double & sigmasqr,
        RngStream & rng);
double RngStream_T(const double df, RngStream rng);
