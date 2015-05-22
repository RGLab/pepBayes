#' @title MCMC sampling for pepBayes models.
#' 
#' @description Sample the posterior distribution of pepBayes models. The generic 
#' function has methods for \link{peptideSet} objects and for \link{pepBayesFit}
#' objects for restarting a Markov chain from its last position.
#' 
#' @details PepBayes estimates the posterior probability of
#' a differential antibody response against each peptide in a peptide microarray
#' assay, for each subject, when compared with control samples. Two common experimental 
#' designs are automatically 
#' detected. In an unpaired design, differences between a population of interest and 
#' a control population are detected. A paired design has matched samples for each subject
#' before and after a treatment is applied. The function \code{pepBayesMcmc} 
#' automatically detects the experimental design from the data, and
#' runs MCMC to sample the posterior distribution of the appropriate pepBayes model.
#' 
#' When \code{x} is a \code{peptideSet} object, the paired or unpaired model is 
#' chosen depending on the 'ptid' column of \code{pData(x)}. Paired data is detected
#' when each unique value of 'ptid' appears exactly twice among all slides.
#' 
#' The \code{position_data} argument supplies information about peptide position
#' information. This argument may be a 'data.frame' with columns 'start',
#' 'end' or 'width', and 'peptide'. If a 'GRanges' object, then it must have either peptide
#' as a name or have peptide as a metadata column. If omitted (\code{position_data = NULL}),
#' response probability hyperparamters will be fixed during estimation.
#' 
#' @param x A pepBayesFit or peptideSet object containing experiment data.
#' @param ... Additional arguments for sampling control. See methods.
#' @return An object of type \link{pepBayesFit} containing sorted slide metadata,
#'  Markov chain posterior samples, posterior probabilities of binding, and
#'  chain information.
#' 
#' @seealso 
#' \code{\link{pepBayesFit}}
#' \code{\link{pepBayesFit-methods}} 
#' \code{\link{peptideSet}} 
#' \code{\link{pepBayesEcm}}
#' \code{\link{pepBayesCalls}}
#' 
#' @docType methods
#' @rdname pepBayesMcmc-methods
#' @exportMethod pepBayesMcmc
#' @useDynLib pepBayes
setGeneric("pepBayesMcmc", function(x, ...) standardGeneric("pepBayesMcmc"))

#' @rdname pepBayesMcmc-methods
#' @param position_data A data.frame or GRanges object with peptide position
#' information.
#' @param control_id Character, indicating which slides are control
#' slides from the "visit" column of the phenoData slot of the peptideSet argument.
#' @param n_samples Number of MCMC posterior samples to collect.
#' @param n_thin Number of MCMC iterations between posterior samples.
#' @param n_burn Number of MCMC iterations to perform before sampling begins.
#' @param n_threads Number of parallel threads to use during sampling.
#' @param prior_control_list A named list of parameter values to initialize
#' @param chain_control_list A named list of parameters controlling 
#' which parameters to update or fix.
#' @aliases pepBayesMcmc,peptideSet-method
setMethod("pepBayesMcmc", signature = "peptideSet",
        function(x, position_data = NULL, control_id = NULL,
                n_samples, n_thin, n_burn, n_threads = detectCores(), 
                prior_control_list = NULL, chain_control_list = NULL) {    
            data_list = .joinPeptideData(x, position_data, control_id)
            par_list = .initializeParameters(data_list, prior_control_list)
            chain_list = .initializeChainList(n_samples, n_thin, n_burn,
                    data_list$method, n_threads, seed = 123, chain_control_list)
            if (is.null(position_data))
                chain_list$update_omega_hypers = FALSE
            message(paste("Detected", data_list$method, "data"))
            result = .Call("pepBayesMcmc", par_list, data_list,
                    chain_list, PACKAGE = "pepBayes")
            result = .processMcmcOutput(result, data_list,
                    par_list, chain_list)
            return(result)
        })

#' @details When \code{x} is a pepBayesFit object, the chain is started from the most
#' recent parameter values in the fitted object. In the event that \code{x} is
#' the output of a previous MCMC fit, the argument \code{join_output} becomes 
#' relevant. If \code{TRUE}, parameter traces are appended to those in \code{x}
#' and the posterior probabilities of binding are averaged together based on
#' the number of iterations in the fitted object and the current chain.
#' @rdname pepBayesMcmc-methods
#' @param join_output A boolean. Should input be merged with output?
#' @aliases pepBayesMcmc,pepBayesFit-method
setMethod("pepBayesMcmc", signature = "pepBayesFit",
        function(x, n_samples, n_thin, n_burn = 0, n_threads = detectCores(),
                join_output = TRUE) {
            if (x@fitType == "ecm") {
                # ecm algorithm is on normal, not logit scale.
                om = x@algState$omega
                om = ifelse(om < 1.0 && om > 0.0, log(om / 1 - om),
                        log((.001 + .99 * om) / (.999 - .99 * om)))
                x@algState$omega = om
            }
            chain_list = .initializeChainList(n_samples, n_thin,
                    n_burn, x@slideInfo$method, n_threads, seed = 123, x@algInfo)
            par_list = .initializeParameters(x@slideInfo, x@algState)
            result = .Call("pepBayesMcmc", par_list, x@slideInfo,
                    chain_list, PACKAGE = "pepBayes")
            result = .processMcmcOutput(result, x@slideInfo, x@algState,
                    chain_list)
            if (join_output && x@fitType == "mcmc")
                result = .joinMcmcOutput(x, result)
            return(result)
        })

#' @title Posterior maximization for pepBayes models.
#' 
#' @description Maximize the posterior distribution of pepBayes models. The generic 
#' function has methods for \link{peptideSet} objects and for \link{pepBayesFit}
#' objects for restarting an optimization from its last position.
#' 
#' @details PepBayes estimates the posterior probability of
#' a differential antibody response against each peptide in a peptide microarray
#' assay, for each subject,
#' when compared with control samples. Two common experimental designs are automatically 
#' detected. In an unpaired design, differences between a population of interest and 
#' a control population are detected. A paired design has matched samples for each subject
#' before and after a treatment is applied. The function \code{pepBayesEcm} 
#' automatically detects the experimental design from the data.
#' 
#' When \code{x} is a \code{peptideSet} object, the paired or unpaired model is 
#' chosen depending on the 'ptid' column of \code{pData(x)}. Paired data is detected
#' when each unique value of 'ptid' appears exactly twice among all slides.
#' 
#' The \code{position_data} argument supplies information about peptide position
#' information. This argument may be a 'data.frame' with columns 'start',
#' 'end' or 'width', and 'peptide'. If a 'GRanges' object, then it must have either peptide
#' as a name or have peptide as a metadata column. Response probability hyperparamters are
#' fixed for EM estimation, but if position data are omitted (\code{position_data = NULL}), 
#' they will also be fixed if the resulting fitted object is used to initialize a 
#' \code{pepBayesMcmc}.
#' 
#' @param x A pepBayesFit or peptideSet object containing experiment data.
#' @param ... Additional arguments for sampling control.
#' @return An object of type \link{pepBayesFit} containing sorted slide metadata,
#'  parameter trace, posterior probabilities of binding, and
#'  iteration information.
#' 
#' @seealso 
#' \code{\link{pepBayesFit}}
#' \code{\link{pepBayesFit-methods}} 
#' \code{\link{peptideSet}}
#' \code{\link{pepBayesMcmc}}
#' \code{\link{pepBayesCalls}}

#'  
#' @docType methods
#' @rdname pepBayesEcm-methods
#' @exportMethod pepBayesEcm
#' @useDynLib pepBayes
setGeneric("pepBayesEcm", function(x, ...) standardGeneric("pepBayesEcm"))

#' @param position_data A data.frame or GRanges object with peptide position
#' information.
#' @param n_iter The number of EM iterations to perform.
#' @param n_threads The number of parallel threads to use in the E-step.
#' @param n_rule The number of nodes to use for Gauss-Hermite numerical integration.
#' @param n_isamp The number of importance samples to estimate the E-step.
#' @param prior_control_list A named list controlling initial parameter values.
#' @param iter_control_list A named list controlling which parameters are updated.
#' @param schedule_iterations A boolean. Should the number of importance samples
#' gradually increase to \code{n_isamp} as iterations progress?
#' @rdname pepBayesEcm-methods
#' @aliases pepBayesEcm,peptideSet-method
setMethod("pepBayesEcm", signature = "peptideSet",
        function(x, position_data = NULL, control_id = NULL,
                n_iter = 30, n_threads = 1,  n_rule = 10, n_isamp = 20, 
                prior_control_list = NULL, iter_control_list = NULL,
                schedule_iterations = TRUE) {
            data_list = .joinPeptideData(x, position_data, control_id,
                    n_rule)
            par_list = .initializeParameters(data_list, prior_control_list)
            iter_list = .initializeIterList(n_iter, n_threads, n_isamp,
                    iter_control_list)
            iter_list$schedule_iterations = schedule_iterations
            if (is.null(position_data))
                iter_list$update_omega_hypers = FALSE
            message(paste("Detected", data_list$method, "data"))
            result = .Call("pepBayesEcm", par_list, data_list, iter_list,
                    PACKAGE = "pepBayes")
            result = .processEcmOutput(result, data_list, par_list,
                    iter_list)
            return(result)
        })

#' @details When \code{x} is a pepBayesFit object, the optimization is started from the most
#' recent parameter values in the fitted object. In the event that \code{x} is
#' the output of a previous ECM fit, the argument \code{join_output} becomes 
#' relevant. If \code{TRUE}, parameter traces are appended to those in \code{x}.
#' @param join_output A boolean. Should input be merged with output?
#' @rdname pepBayesEcm-methods
#' @aliases pepBayesEcm,pepBayesFit-method
setMethod("pepBayesEcm", signature = "pepBayesFit",
        function(x, n_iter, n_threads = 1, n_rule = 10, n_isamp = 20,
                schedule_iterations = FALSE, join_output = TRUE) {
            if (x@fitType == "mcmc") {
                # handle transformation of omega between algorithms.
                om = 1 / (1 + exp(-x@algState$omega))
                x@algState$omega = om                                
            }
            iter_list = .initializeIterList(n_iter, n_threads, n_isamp,
                    x@algInfo)
            iter_list$schedule_iterations = schedule_iterations
            par_list = .initializeParameters(x@slideInfo, x@algState)
            result = .Call("pepBayesEcm", par_list, x@slideInfo,
                    iter_list, PACKAGE = "pepBayes")
            result = .processEcmOutput(result, x@slideInfo, x@algState,
                    iter_list)
            if (join_output && x@fitType == "ecm") {
                result = .joinEcmOutput(x, result)
            }
            return(result)
        })

#' @title pepBayesFit methods.
#' 
#' @description Methods for interacting with \code{\link{pepBayesFit}} objects.
#' 
#' @section Accessors:
#' \describe{
#'  \item{\code{position(x)}}{Return the position of peptides.}
#'  \item{\code{peptide(x)}}{Return the peptide amino acid sequences.}
#'  \item{\code{parTrace(x)}}{Return trace of fitted parameter values across iterations
#' of fitting algorithm.}
#'  \item{\code{ppb(x)}}{Return the estimated posterior probability of binding for each subject
#' /peptide combination.}
#' }
#' 
#' @section pepBayes Analysis:
#' \describe{
#'  \item{\code{pepBayesCalls(x, fdr)}}{Return calls for a pepBayesFit at the desired
#' \code{fdr} value.}
#'  \item{\code{pepBayesMcmc(x, ...)}}{Run \code{\link{pepBayesMcmc}}, starting the chain
#' at the last position of the parameters in \code{x}.}
#'  \item{\code{pepBayesEcm(x, ...)}}{Run \code{\link{pepBayesEcm}}, starting the optimization
#' at the last position of the parameters in \code{x}.}
#' }
#' 
#' @seealso 
#' \code{\link{pepBayesEcm}}
#' \code{\link{pepBayesMcmc}}
#' \code{\link{pepBayesCalls}}
#' 
#' @name pepBayesFit-methods
#' @rdname pepBayesFit-methods
#' @exportMethod "position"
#' @exportMethod "peptide"
#' @exportMethod "parTrace"
#' @exportMethod "ppb"
NULL

setMethod("position", signature = "pepBayesFit",
        function(x) {
            return(x@slideInfo$orig_pos)
        })

setMethod("peptide", signature = "pepBayesFit",
        function(x) {
            return(x@slideInfo$peptide)
        })

setGeneric("parTrace", function(x, ...) standardGeneric("parTrace"))

setMethod("parTrace", signature = "pepBayesFit",
        function(x) {
            return(x@parameterOutput)
        })

setGeneric("ppb", function(x, ...) standardGeneric("ppb"))
setMethod("ppb", signature = "pepBayesFit",
        function(x) {
            return(x@ppb)
        })

#' @title Call peptide binding from pepBayes output.
#' 
#' @description Use posterior probability of response to generate a 
#' call threshold that controls the false discovery rate of pepBayes positivity
#' calls.
#' 
#' @details Probabilities are sorted in descending order. The largest integer K
#' is selected such that the mean of the first K sorted probability values exceeds
#' the value (1 - \code{fdr}). The probability threshold is chosen to be the K-th
#' sorted probability value. Entries with probability greater than or equal to the
#' threshold value are considered positive calls, else an entry is a negative call.
#' 
#' @param x A matrix or pepBayesFit object.
#' @param fdr Desired false discovery rate.
#' 
#' @return A matrix of booleans with TRUE denoting a positive call, and FALSE 
#' a negative call.
#' 
#' @rdname pepBayesCalls-methods
#' @docType methods
#' @exportMethod pepBayesCalls
#' @seealso 
#' \code{\link{pepBayesEcm}}
#' \code{\link{pepBayesMcmc}}
#' \code{\link{pepBayesFit-methods}}
setGeneric("pepBayesCalls", function(x, fdr, ...) standardGeneric("pepBayesCalls"))

#' @rdname pepBayesCalls-methods
#' @aliases pepBayesCalls,pepBayesFit,numeric-method
setMethod("pepBayesCalls", signature = c("pepBayesFit", "numeric"),
        function(x, fdr = .05) {
            return(.getCalls(ppb(x), fdr))
        })

#' @rdname pepBayesCalls-methods
#' @aliases pepBayesCalls,matrix,numeric-method
setMethod("pepBayesCalls", signature = c("matrix", "numeric"),
        function(x, fdr = .05) {
            return(.getCalls(x, fdr))
        })





