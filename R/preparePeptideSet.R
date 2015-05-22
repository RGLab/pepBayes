#### functions to prepare a peptideSet object for ECM and MCMC computations
#' @title Check peptideSet phenoData.
#' 
#' @description Check peptideSet phenoData slot for appropriate column names.
#' 
#' @details Internal function.
#' @param peptide_set a peptideSet object.
#' @return NULL
#' 
#' @author Gregory Imholte
#' @keywords internal
.checkPeptideSet = function(peptide_set) {
    pdata_names = names(pData(peptide_set))
    if (!all(c("ptid", "visit") %in% pdata_names)) {
        stop("the phenoData slot of the peptide_set object must contain\n columns \"ptid\" and \"visit\"")
    }
    return(NULL)
}

#' @title Check validity of control_id.
#' 
#' @description Check whether control_id is a valid string found in the
#' phenoData slot of the peptideSet object.
#' 
#' @details Internal function.
#' @param peptide_set a peptideSet object.
#' @param control_id character string indicating which slides are control slides.
#' @return NULL
#' 
#' @author Gregory Imholte
#' @keywords internal
.checkControlId = function(peptide_set, control_id) {
    if (length(control_id) != 1) {
        stop("invalid control_id") 
    }
    if (!(control_id %in% pData(peptide_set)$visit)) {
        stop("control_id not found in visit column of phenoData slot from peptide_set object")
    }
}

#' @title Detect experimental design
#' 
#' @description Detect either the paired or unpaired design for pepBayes modeling.
#' 
#' @details Internal function.
#' @param peptide_set a peptideSet object
#' @return either "paired" or "unpaired"
#' 
#' @author Gregory Imholte
#' @keywords internal
.detectMethod = function(peptide_set) {
    ptid = pData(peptide_set)$ptid
    if (all(table(ptid) == 2)) {
        return("paired")
    }
    return("unpaired")
}


#' @title Reorder peptideSet columns.
#' 
#' @description Reorder peptideSet columns so that data matches expected input
#' to pepBayes algorithms.
#' 
#' @details Internal function. pepBayes expects paired data to be sorted first
#' by ptid, then by pre/post. Unpaired data is sorted first by controls, then by 
#' ptid.
#' @param peptide_set a peptideSet object .
#' @param control_id string to identify control samples.
#' @param method character. "paired" or "unpaired".
#' @return the peptide_set with its columns reordered. 
#' 
#' @author gimholte
#' @keywords internal
.orderPeptideSetColumns = function(peptide_set, control_id, method) {
    pre = control_id
    post = setdiff(unique(pData(peptide_set)$visit), pre)
    visitnum = match(pData(peptide_set)$visit, c(pre, post))
    if (method == "paired") {
        o = with(pData(peptide_set), order(ptid, visitnum))
        return(peptide_set[,o])
    }
    
    o = with(pData(peptide_set), order(visitnum, ptid))
    return(peptide_set[,o])
}

# compute model summary statistics 
#' @title Generate summary statistics for peptideSet.
#' 
#' @description compute summary statistics for peptideSet object, such as mean and 
#' sum of squares.
#' 
#' @details Internal function. Missing values are removed.
#' @param peptide_set 
#' @return A list with peptide, mean, sum of squares, and number or replicates. 
#' 
#' @author Gregory Imholte
#' @keywords internal
.getSummaryStatistics = function(peptide_set) {
    # convert expression to data table for easy processing
    expr_dat = data.table(exprs(peptide_set))
    expr_dat$peptide = peptide(peptide_set)
    setkey(expr_dat, peptide)
    # apply mean and SS function to all columns, by peptide
    y_mean = expr_dat[, lapply(.SD, mean, na.rm = TRUE), by = "peptide"]
    y_ss = expr_dat[, lapply(.SD, function(x) 
                        sum((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE)),
            by = "peptide"]
    # use expression magic to get number of replicates for each peptide
    # note the use of ` to handle spaces in the column names of the
    # expressions.
    expr = parse(text = paste0("list(n_rep = length(`",
                            colnames(expr_dat)[1] ,
                            "`))"))
    y_n_rep = expr_dat[, eval(expr), by = "peptide"]

    # because peptide was the table key, we are assured that
    # y_ss, y_mean, and y_n_rep share a common ordering by
    # peptide
    pep_seq = y_mean[, peptide]
    y_mean[, peptide:= NULL]
    y_n_rep[, peptide:= NULL]
    y_ss[, peptide:= NULL]

  
    # use the data.matrix function to efficiently convert the data
    # tables to numeric matrices
    summary_stats = list(peptide = pep_seq, y_mean = data.matrix(y_mean),
            y_ss = data.matrix(y_ss),
            n_rep = as.double(y_n_rep$n_rep))
    return(summary_stats)
}

# function that creates a lookup table for peptides and their positions,
# sorted by position

#' @title Process position data.
#' 
#' @description Organize position data into a common format understandable
#' to pepBayes algorithms.
#' 
#' @details Internal function. Relies on pepStat function 
#' \code{\link[pepStat]{create_db}}.
#' @param position_data a GRanges or data.frame containing position info 
#' @return Data table with elements peptide and position, ordered by position.
#' 
#' @author Gregory Imholte
.processPositionData = function(position_data, peptide_set) {
    if (is.null(position_data)) {
        return(.defaultPositionInfo(peptide_set))
    }
    pos_db = pepStat::create_db(position_data)
    pos = round((start(ranges(pos_db)) + end(ranges(pos_db))) / 2)
    o = order(pos)
    pep = names(pos_db)
    return(data.table(peptide = pep[o], position = pos[o]))
}

#' @title Create default position data.
#' 
#' @description Create default position data when the user declines to provide it.
#'
#' @details Internal function. Each peptide is assigned a unique position value. 
#' @param peptide_set a peptideSet object.
#' @return  A data.table with elements peptide and position, ordered by position.
#' 
#' @author Gregory Imholte
#' @keywords internal
#' @import data.table
.defaultPositionInfo = function(peptide_set) {
    peps = unique(peptide(peptide_set))
    npep = length(peps)
    return(data.table(peptide = peps, position = 1:npep))
}

#' @title Join peptideSet and position data for pepBayes.
#' 
#' @description Compute summary statitics, order data, join peptide information
#' to data, and create indexing vectors for pepBayes algorithms.
#' 
#' @details Internal function. 
#' 
#' @param peptide_set a peptideSet object 
#' @param position_data position data for peptideSet objects
#' @param control_id character identifying control samples
#' @param n_rule number of gauss hermite quadrature nodes
#' @return A list containing matched summary statistics and position data to be fed
#' into pepBayes algorithms.
#' 
#' @author Gregory Imholte
#' @keywords internal
#' @import data.table
#' @importFrom fastGHQuad gaussHermiteData
.joinPeptideData = function(peptide_set, position_data = NULL, control_id,
        n_rule = 1) {
    # order the columns
    .checkPeptideSet(peptide_set)
    .checkControlId(peptide_set, control_id)
    method = .detectMethod(peptide_set)
    peptide_set = .orderPeptideSetColumns(peptide_set, control_id, method)
    pset_summary = .getSummaryStatistics(peptide_set)
    pos_table = .processPositionData(position_data, peptide_set)
    if (is.null(position_data))
        pset_summary$pos_missing = TRUE
    else
        pset_summary$pos_missing = FALSE
    
    # "not missing" index
    nm_idx = pset_summary$peptide %in% pos_table[, peptide]
    if (!all(nm_idx)) {
        n_missing = sum(!(nm_idx))
        warning(paste(n_missing, 
                        "peptides in peptide_set were not found in position_data and are\nomitted from analysis"))
    }
    lnames = names(pset_summary)
    pset_summary = lapply(pset_summary, function(x, idx) {
                d = dim(x)
                x = x[idx]
                if (!is.null(d)) {
                    dim(x) = c(sum(idx), d[2])
                }
                return(x)
            }, idx = nm_idx)
    names(pset_summary) = lnames
    
    # reorder peptides and summarized data by position
    match_idx = match(pos_table[, peptide], pset_summary$peptide)
    pset_summary = lapply(pset_summary, function(x, idx) {
                d = dim(x)
                if (is.null(d)) {
                    return(x[idx])
                }
                return(x[idx,])
            }, idx = match_idx)
    names(pset_summary) = lnames
    
    unique_pos = unique(pos_table[, position])
    pos_original_scale = pos_table[, position][match(pset_summary$peptide, pos_table[, peptide])]
    pos = rep(1:length(unique_pos), times = table(pos_table[, position])) - 1
    pos_count = as.vector(table(pos_table[, position]))
    n_pos = length(unique_pos)
    pos_start = c(0, (cumsum(pos_count) - 1)[1:(n_pos - 1)])
    
    pset_summary$pos = pos
    pset_summary$orig_pos = pos_original_scale
    pset_summary$n_position = as.integer(n_pos)
    pset_summary$pos_count = pos_count
    pset_summary$pos_start = pos_start
    pset_summary$method = method
    pset_summary$metadata = pData(peptide_set)
    
    ctl_ind = pData(peptide_set)$visit == control_id
    pset_summary$ctl_ind = ctl_ind
    if (method == "unpaired") {
        # divide trt and ctl samples
        pset_summary$y_mean_trt = pset_summary$y_mean[, !ctl_ind]
        pset_summary$y_mean_ctl = pset_summary$y_mean[, ctl_ind]
        pset_summary$y_mean = NULL
        
        pset_summary$y_ss_trt = pset_summary$y_ss[, !ctl_ind]
        pset_summary$y_ss_ctl = pset_summary$y_ss[, ctl_ind]
        pset_summary$y_ss = NULL
    }
    
    gauss_rule = fastGHQuad::gaussHermiteData(n_rule)
    pset_summary$gauss_rule_x = gauss_rule$x
    pset_summary$gauss_rule_w = log(gauss_rule$w)
    return(pset_summary)
}















