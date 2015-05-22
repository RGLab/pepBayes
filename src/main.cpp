/*
 * main.cpp
 *
 *  Created on: Oct 22, 2013
 *      Author: hoblitz
 */

#include "RcppArmadillo.h"
#include "RngStream.h"
#include "Rv.h"
#include <omp.h>
#include "Densities.h"
#include "AdaptiveGibbsSampler.h"
#include "FitBase.h"
#include "PairedMcmc.h"
#include "PairedEcm.h"
#include "UnpairedMcmc.h"
#include "UnpairedEcm.h"

#include <string.h>

RcppExport SEXP pepBayesMcmc(SEXP r_par_list, SEXP r_data_list, SEXP r_chain_list) {
    Rcpp::List par_list(r_par_list);
    Rcpp::List data_list(r_data_list);
    Rcpp::List chain_list(r_chain_list);
    if (Rcpp::as<string>(data_list["method"]) == "unpaired") {
        UnpairedMcmc chain(par_list, data_list, chain_list);
        return chain.run();
    } else {
        PairedMcmc chain(r_par_list, r_data_list, r_chain_list);
        return chain.run();
    }
}

RcppExport SEXP pepBayesEcm(SEXP r_par_list,
        SEXP r_data_list, SEXP r_iter_list) {
    Rcpp::List par_list(r_par_list);
    Rcpp::List data_list(r_data_list);
    Rcpp::List iter_list(r_iter_list);
    if (Rcpp::as<string>(data_list["method"]) == "unpaired") {
        UnpairedEcm ecm(par_list, data_list, iter_list);
        return ecm.run();
    } else {
        PairedEcm ecm(par_list, data_list, iter_list);
        return ecm.run();
    }
}

