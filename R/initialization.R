#' @title Default prior values for pepBayes models.
#' 
#' @description Return a list with default prior for pepBayes models.
#' 
#' @details Internal function
#' 
#' @return A named list with default values for pepBayes model hyperparameters.
#' @keywords internal
.defaultPrior = function() {
    return(list(
        "m_beta0_prior_mean" = 4.0,
        "m_beta0_prior_prec" = .25,
        "m_beta1_prior_mean" = 2.0,
        "m_beta1_prior_prec" = 2.0,
        "nu_beta0_shape" = .5,
        "nu_beta0_rate" = .5,
        "nu_beta1_shape" = .5,
        "nu_beta1_rate" = .5,
        "lambda_alpha_rate" = 1000,
        "lambda_eps_rate" = 1000,
        "s_alpha_rate" = 1000,
        "s_eps_rate" = 1000,
        "nu" = 4,
        "lambda_a" = 1.0,
        "lambda_b" = .1))
}

#' @title Check replacement prior value.
#' 
#' @description Check validity of parameter values to replace.
#' 
#' @details Internal function
#' 
#' @param character, name of the parameter
#' @param cur_par_val, current value of parameter
#' @param new_par_val, value to replace current value
#' 
#' @return A boolean indicating whether the new parameter value
#' is of the same length. In the event of failure, the "message" attribute
#' of the return value has a reason for failure.
#' @keywords internal
.checkNewPriorValue = function(par_name, cur_par_val, new_par_val) {
    is_ok = TRUE
    e_message = NULL
    if (!is.numeric(new_par_val)) {
        e_message = "new parameter values must be numeric"
        is_ok = FALSE
    }
    if (length(cur_par_val) != length(new_par_val)) {
        e_message = "new parameter length must match current parameter length"
        is_ok = FALSE
    }
    attr(is_ok, "message") = e_message
    return(is_ok)
}

#' @title Update the prior with new values.
#' 
#' @description Update the parameter list based on a list of new values.
#' 
#' @details Internal function
#' @param cur_par_list named list of current parameter values
#' @param new_prior_values named list of parameters to replace, with values
#' 
#' @return The list of parameter values, with values in new_prior_values
#' replacing corresponding entries in cur_par_list.
#' @keywords internal
.updateDefaultPrior = function(cur_par_list, new_prior_values = NULL) {
    if (is.null(new_prior_values))
        return(cur_par_list)
    
    if (!is.list(new_prior_values))
        new_prior_values = as.list(new_prior_values)
    
    npv_names = names(new_prior_values)
    cpl_names = names(cur_par_list)
    
    if (is.null(npv_names))
        stop("new values for prior must be a named vector or list")
    
    if (!all(npv_names %in% cpl_names)) {
        bad_ind = which(!(npv_names %in% cpl_names))
        bad_names = npv_names[bad_ind]  
        e_message = paste(paste(bad_names, collapse = ", "),
                "parameter(s) do not match prior parameter names")
        stop(e_message)
    }
    
    for (par_name in npv_names) {
        new_par_val = new_prior_values[[par_name]]
        
        is_ok = .checkNewPriorValue(par_name, 
                cur_par_list[[par_name]], new_par_val)
        if (!is_ok) stop(attr(is_ok, "message"))
        
        cur_par_list[[par_name]] = new_par_val
    }
    return(cur_par_list)
}

#' @title Initialize pepBayes parameter values.
#' 
#' @description Initialize pepBayes parameter values based on data.
#' 
#' @details Internal function. 
#' 
#' @param data_list result of .joinPeptideData, contains all probe summaries,
#' slide info, and position info.
#' @param prior_control_list named list of parameter values to set.
#' 
#' @return A named list with initial parameter values for starting MCMC/ECM.
#' @keywords internal
.initializeParameters = function(data_list, prior_control_list) {
    if (data_list$method == "paired") {
        pre_ind = seq(1, ncol(data_list$y_mean) - 1, by = 2)
        n_slide = ncol(data_list$y_mean)
        mu_init = colMeans(data_list$y_mean)
        y_mean = data_list$y_mean[, pre_ind]
        y_ss = data_list$y_ss[, pre_ind]
    } else {
        y_mean = data_list$y_mean_ctl
        y_ss = data_list$y_ss_ctl
        mu_init = with(data_list, colMeans(cbind(y_mean_ctl, y_mean_trt)))
        n_slide = with(data_list, ncol(y_mean_trt) + ncol(y_mean_ctl))
    }

    n_peptide = nrow(y_mean)
    n_position = data_list$n_position
    par_list = c(.defaultPrior(), list(beta0 = rowMeans(y_mean),
                    beta1 = rep(2.0, n_peptide),
                    nu_eps = (data_list$n_rep * (ncol(y_ss) - 1)) / 
                            rowSums(y_ss),
                    nu_alpha0 = 1 / apply(y_mean, 1, var),
                    nu_alpha1 = 1 / apply(y_mean, 1, var),
                    mu = mu_init - mean(mu_init),
                    omega = rep(.05, n_peptide),
                    a = rep(.125, n_position),
                    b = rep(1.125, n_position),
                    m_beta0 = 2.0,
                    m_beta1 = 2.0,
                    nu_beta0 = 2.0,
                    nu_beta1 = 2.0,
                    lambda_eps = .05,
                    s_eps = .05,
                    lambda_alpha = .05,
                    s_alpha = .05))
    
    if (data_list$method == "paired") {
        par_list$gamma = matrix(0.0, n_peptide, n_slide / 2)
        par_list$alpha0 = matrix(0.0, n_peptide, n_slide / 2)
        par_list$alpha1 = matrix(0.001, n_peptide, n_slide / 2) * par_list$gamma
        par_list$w = matrix(1.0, n_peptide, n_slide)
        par_list$u = matrix(1.0, n_peptide, n_slide)
    } else {
        n_trt = ncol(data_list$y_mean_trt)
        n_ctl = ncol(data_list$y_mean_ctl)
        par_list$alpha_trt = matrix(0.0, n_peptide, n_trt)
        par_list$alpha_ctl = matrix(0.0, n_peptide, n_ctl)
        par_list$gamma = matrix(0.0, n_peptide, n_trt)
        par_list$w_trt = matrix(1.0, n_peptide, n_trt)
        par_list$w_ctl = matrix(1.0, n_peptide, n_ctl)
        par_list$u_trt = matrix(1.0, n_peptide, n_trt)
        par_list$u_ctl = matrix(1.0, n_peptide, n_ctl)
    }
    
    S = diag(n_slide) - 1 / n_slide
    par_list$Q = eigen(S, symmetric = TRUE)$vector[,1:(n_slide - 1)]
    par_list = .updateDefaultPrior(par_list, prior_control_list)
    return(par_list)
}

#' @title Initialize MCMC control.
#' 
#' @description Initialize parameters controling MCMC updates and sample size.
#' 
#' @param n_samples number of posterior samples to collect.
#' @param n_thin number of iterations between collected samples
#' @param n_burn number of burn-in iterations
#' @param method character, equal to "paired" or "unpaired"
#' @param n_threads number of parallel processing threads to use
#' @param seed seed for pepBayes random number generator
#' @param chain_pars_update named list of chain parameters to change from
#' default values
#' 
#' @details Internal function
#' 
#' @return A list of values indicating sampling scheme, and which
#' model parameters to fix or update
#' @keywords internal
.initializeChainList = function(n_samples, n_thin, n_burn, method,
        n_threads, seed, chain_pars_update) {
    chain_list = list("n_samples" = n_samples, "n_thin" = n_thin,
            "n_burn" = n_burn, "method" = method, "n_threads" = n_threads,
            "seed" = seed)
    update_list = list(
            "update_mu" = TRUE,
            "update_beta0" = TRUE,
            "update_beta1" = TRUE,
            "update_beta_hypers" = TRUE,
            "update_precision" = TRUE,
            "update_precision_hypers" = TRUE,
            "update_omega" = TRUE,
            "update_omega_hypers" = TRUE,
            "update_alpha" = TRUE,
            "update_weights" = TRUE,
            "update_gamma" = TRUE)
    
    if (is.null(chain_pars_update)) {
        return(c(chain_list, update_list))
    }
    cur_names = names(update_list)
    new_names = names(chain_pars_update)
    changed = intersect(cur_names, new_names)
    update_list[changed] = chain_pars_update[changed]
    chain_list = c(chain_list, update_list)
    return(chain_list)
}

#' @title Initialize EM control.
#' 
#' @description Initialize control for EM algorithm, such as which model
#' parameters to update and how many iterations to perform.
#' 
#' @details Internal function.
#' @param n_iter number of EM iterations to perform.
#' @param n_threads number of parallel threads to use.
#' @param n_isamp number of importance samples to employ.
#' @param iter_pars_update named list of new chain control values to update. 
#' @return A named list of algorithm control values.
#' 
#' @author Gregory Imholte
#' @keywords internal
.initializeIterList = function(n_iter,
        n_threads, n_isamp = 10, iter_pars_update) {
    iter_list = list("n_iter" = n_iter, "n_threads" = n_threads,
            "n_isamp" = n_isamp)
    update_list = list(
            "update_mu" = TRUE,
            "update_beta0" = TRUE,
            "update_beta1" = TRUE,
            "update_beta_hypers" = TRUE,
            "update_precision" = TRUE,
            "update_precision_hypers" = TRUE,
            "update_omega" = TRUE,
            "update_omega_hypers" = TRUE,
            "update_weights" = TRUE)
    if (is.null(iter_pars_update)) {
        return(c(iter_list, update_list))
    }
    cur_names = names(update_list)
    new_names = names(iter_pars_update)
    changed = intersect(cur_names, new_names)
    update_list[changed] = iter_pars_update[changed]
    iter_list = c(iter_list, update_list)
    return(iter_list)
}

