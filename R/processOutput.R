#' @title Create a pepBayesFit object from MCMC output.
#' 
#' @description Combine user inputs and MCMC output into a pepBayesFit object containing all
#' model output, as well as data summaries and user inputs required to restart the chain.
#' 
#' @details Internal function.
#' 
#' @param output_list output from pepBayes C++ engine.
#' @param data_list output from .joinPeptideData containing ordered metadata,
#'  summarized data, position information, gauss hermite quadrature info.
#' @param par_list output from .initializeParameter, containing user updates to 
#'  prior parameter values and mu eigendecomposition.
#' @param chain_list output from .initializeChainList setting more obscure arguments
#'  that control which model parameters are updated.
#' @return A pepBayesFit object with all model output.
.processMcmcOutput = function(output_list,
        data_list, par_list, chain_list) {
    # give parameter traces appropriate dimnames
    rownames(output_list$ppb) = data_list$peptide
    rownames(output_list$mu) = with(data_list, 
            paste(metadata$ptid, metadata$visit, sep = "_"))
    rownames(output_list$beta0) = data_list$peptide
    rownames(output_list$beta1) = data_list$peptide
    rownames(output_list$nu_eps) = data_list$peptide
    rownames(output_list$nu_alpha0) = data_list$peptide
    rownames(output_list$nu_alpha1) = data_list$peptide
    rownames(output_list$omega) = data_list$peptide
    rownames(output_list$hypers) = c("s_alpha", "lambda_alpha",
            "s_eps", "lambda_eps", "m_beta0", "nu_beta0", 
            "m_beta1", "nu_beta1", "nu_err", "nu_re")
    colnames(output_list$ppb) = with(data_list, metadata$ptid[!ctl_ind])
    
    par_output = output_list[c("hypers", "beta0", "beta1", "mu", "nu_eps",
                    "nu_alpha0", "nu_alpha1", "omega", "a", "b")]
    ppb = output_list[["ppb"]]
    
    # extract chain state at last iteration
    n_samples = chain_list$n_samples
    chain_state = with(output_list,
            list(beta0 = beta0[, n_samples],
                    beta1 = beta1[, n_samples],
                    mu = mu[, n_samples],
                    nu_eps = nu_eps[, n_samples],
                    nu_alpha0 = nu_alpha0[, n_samples],
                    nu_alpha1 = nu_alpha1[, n_samples],
                    omega = omega[, n_samples],
                    a = a[, n_samples],
                    b = b[, n_samples],
                    gamma = gamma))
    chain_state[names(.defaultPrior())] = par_list[names(.defaultPrior())]
    chain_state[rownames(output_list$hypers)] = 
            as.list(output_list$hypers[, n_samples])
    chain_state[["Q"]] = par_list$Q
    if (data_list$method == "unpaired") {
        varnames = c("u_trt", "u_ctl", "w_trt", "w_ctl",
                "alpha_trt", "alpha_ctl")
    } else {
        varnames = c("u", "w", "alpha0", "alpha1")
    }
    chain_state[varnames] = output_list[varnames]
    
    pFit = new("pepBayesFit", parameterOutput = par_output, ppb = ppb,
            slideInfo = data_list, algInfo = chain_list,
            fitType = "mcmc", algState = chain_state)
    return(pFit)
}

#' @title Create a pepBayesFit object from ECM output.
#' 
#' @description Combine user inputs and ECM output into a pepBayesFit object containing all
#' model output, as well as data summaries and user inputs required to 
#' restart the algorithm from it's most recent iteration.
#' 
#' @details Internal function.
#' 
#' @param output_list output from pepBayes C++ engine.
#' @param data_list output from .joinPeptideData containing ordered metadata,
#'  summarized data, position information, gauss hermite quadrature info.
#' @param par_list output from .initializeParameter, containing user updates to 
#'  prior parameter values and mu eigendecomposition.
#' @param iter_list output from .initializeIterList setting more obscure arguments
#'  that control which model parameters are updated.
#' @return A pepBayesFit object with all model output.
.processEcmOutput = function(output_list,
        data_list, par_list, iter_list) {
    # give parameter traces appropriate dimnames
    rownames(output_list$ppb) = data_list$peptide
    rownames(output_list$mu) = with(data_list, 
            paste(metadata$ptid, metadata$visit, sep = "_"))
    rownames(output_list$beta0) = data_list$peptide
    rownames(output_list$beta1) = data_list$peptide
    rownames(output_list$nu_eps) = data_list$peptide
    rownames(output_list$nu_alpha0) = data_list$peptide
    rownames(output_list$nu_alpha1) = data_list$peptide
    rownames(output_list$omega) = data_list$peptide
    rownames(output_list$hypers) = c("s_alpha", "lambda_alpha",
            "s_eps", "lambda_eps", "m_beta0", "nu_beta0", 
            "m_beta1", "nu_beta1")
    colnames(output_list$ppb) = with(data_list, metadata$ptid[!ctl_ind])
    
    par_output = output_list[c("hypers", "beta0", "beta1", "mu", "nu_eps",
                    "nu_alpha0", "nu_alpha1", "omega")]
    ppb = output_list[["ppb"]]
    par_info = par_list
    
    # extract chain state at last iteration
    n_iter = iter_list$n_iter
    chain_state = with(output_list,
            list(beta0 = beta0[, n_iter],
                    beta1 = beta1[, n_iter],
                    mu = mu[, n_iter],
                    nu_eps = nu_eps[, n_iter],
                    nu_alpha0 = nu_alpha0[, n_iter],
                    nu_alpha1 = nu_alpha1[, n_iter],
                    omega = omega[, n_iter]))
    chain_state[names(.defaultPrior())] = par_list[names(.defaultPrior())]
    chain_state[rownames(output_list$hypers)] = 
            as.list(output_list$hypers[, n_iter])
    chain_state[c("a", "b", "Q")] = par_list[c("a", "b", "Q")]

    pFit = new("pepBayesFit", parameterOutput = par_output, ppb = ppb,
            slideInfo = data_list, algInfo = iter_list,
            fitType = "ecm", algState = chain_state)
    return(pFit)
}

#' @title Join pepBayes fits from ECM output.
#' 
#' @description Concatenate the parameter trace from ECM outputs.
#' 
#' @details Internal function.
#' @param old_fit pepBayesFit object used to restart algorithm.
#' @param new_fit pepBayesFit output from a restart.
#' @return a pepBayesFit object with joined parameter traces. 
#' 
#' @author Gregory Imholte
#' @keywords internal
.joinEcmOutput = function(old_fit, new_fit) {
    output_names = c("mu", "beta0", "beta1", "nu_eps", "nu_alpha0",
            "nu_alpha1", "omega", "hypers")
    for(n in output_names) {
        new_fit@parameterOutput[[n]] = cbind(old_fit@parameterOutput[[n]],
                new_fit@parameterOutput[[n]])
    }
    
    new_fit@algInfo$n_iter = old_fit@algInfo$n_iter + new_fit@algInfo$n_iter
    return(new_fit)
}

#' @title Join pepBayes fits from MCMC output.
#' 
#' @description Concatenate the parameter trace from MCMC outputs and
#' update the posterior probability of binding.
#' 
#' @details Internal function.
#' @param old_fit pepBayesFit object used to restart algorithm.
#' @param new_fit pepBayesFit output from a restart.
#' @return a pepBayesFit object with joined parameter traces. 
#' 
#' @author Gregory Imholte
#' @keywords internal
.joinMcmcOutput = function(old_fit, new_fit) {
    output_names = c("mu", "beta0", "beta1", "nu_eps", "nu_alpha0",
            "nu_alpha1", "omega", "hypers", "a", "b")
    for(n in output_names) {
        new_fit@parameterOutput[[n]] = cbind(old_fit@parameterOutput[[n]],
                new_fit@parameterOutput[[n]])
    }
    n_samp_old = old_fit@algInfo$n_samples
    n_samp_new = new_fit@algInfo$n_samples
    n_samp_total = n_samp_old + n_samp_new
    
    new_fit@ppb = ppb(new_fit) * (n_samp_new / n_samp_total) +
            ppb(old_fit) * (n_samp_old / n_samp_total)
    
    new_fit@algInfo$n_samples = n_samp_total
    return(new_fit)
}

#' @title Implementation of pepBayesCalls calling method.
#' 
#' @description Estimate whether a peptide experienced true binding.
#' 
#' @details Internal function.
#' @param ppb_mat matrix of posterior probability of binding.
#' @param fdr numeric, false discovery rate threshold.
#' @return a boolean matrix indicating which peptide/subject combinations
#' are called positive/negative. 
#' 
#' @author Gregory Imholte
#' @keywords internal
.getCalls = function(ppb_mat, fdr = .05) {
    o = order(ppb_mat, decreasing = TRUE)
    idx = which(cumsum(ppb_mat[o]) / 1:length(ppb_mat) >= 1 - fdr)
    calls = matrix(FALSE, nrow = nrow(ppb_mat), ncol = ncol(ppb_mat))
    dimnames(calls) = dimnames(ppb_mat)
    if (length(idx) == 0)
        return(calls)
    calls_idx = o[1:max(idx)]
    calls[calls_idx] = TRUE
    attr(calls, "threshold") = ppb_mat[o[max(idx)]]
    return(calls)
}
