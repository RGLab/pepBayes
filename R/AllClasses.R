#' Class pepBayesFit.
#' 
#' Class \code{pepBayesFit} holds the output of pepBayes models and
#' the sorted, summarized data. A handful of accessor functions are
#' implemented to access object slots. See \code{\link{pepBayesFit-methods}}.
#' 
#' @section Slots:
#' \describe{
#'  \item{parameterOutput}{A named list with trace of parameter values
#'  generated as the fit progressed.}
#'  \item{ppb}{A numeric matrix containing estimates of the posterior
#' probabilities of binding for eac subject and peptide combination.}
#'  \item{fitType}{Character, "mcmc" or "ecm".}
#'  \item{slideInfo}{A named list whose entries are various
#' data summaries and internal quantities required for the fit.}
#'  \item{algInfo}{A named list whose values indicate which parameters
#' were fixed, how many iterations were performed, etc.}
#'  \item{algState}{A list with parameter values at the most
#' recent iteration. Used to initialize or restart another pepBayes fit.}
#' }
#' @name pepBayesFit
#' @rdname pepBayesFit-class
#' @exportClass pepBayesFit
#' @aliases pepBayesFit-class
setClass("pepBayesFit",
        representation(parameterOutput = "list", ppb = "matrix",
                fitType = "character", slideInfo = "list", algInfo = "list", 
                algState = "list"))
