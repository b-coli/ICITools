#' Compute specification table for cell-type expression data
#'
#' @description This function implements the Index of Cell Identity speficiation
#' table computation algorithm developed by Idan Efroni and colleagues (Efroni,
#' et al., 2015). This package takes advantage of the \code{furrr::future_map}
#' utility, which enables parallelization on multicore machines. To enable this
#' functionality, you must specify,
#' \code{future::plan(strategy = "multisession")} command prior to executing
#' ICITools functions.

#' @param expression_data a data frame (or tibble) containing the following
#' columns for each gene/dataset combination:
#'
#' \itemize{
#'   \item Locus (the gene or probeset)
#'   \item Expression (the normalized expression value, not log-transformed)
#'   \item Cell_Type (the cell type from which the expression data came from)
#'   \item Sample_Name (the sample where Locus/Cell_Type Expression value was
#'   measured.)
#' }
#'
#' It is assumed that the expression_data object contains no missing values.
#' This is important, since the specificity score computation should be
#' comparable between loci for the same cell types, which would not be the case
#' if some loci/cell-type combinations are missing. As an initial pre-processing
#' step, this function will remove any loci that have missing values.
#'
#' @param bin_method character (only implemented method is "Efroni") or
#' user-defined method for binning expression data. The user-defined method must
#' take in a data.frame with columns, "Cell_Type" and "Expression" for a single
#' locus, and return a data.frame containing columns "Cell_Type", "Expression",
#' and "bin". If \code{bin_method} is set to something other than a function or
#' "Efroni", this function will exit with an error.
#'
#' @param mean_method character (only implemented methods are "Efroni" and "median")
#' or user-defined method for computing expression mean for each locus/cell
#' type. The user-defined method must take in a \code{data.frame} with columns,
#' "Cell_Type" and "Expression", and return a \code{data.frame} containing
#' columns "Cell_Type", "Expression", and "mean_expr". There should be one mean
#' value returned for each Cell_Type/Locus combination. If mean_method is set to
#' something other than a function or "Efroni", mean expression for each
#' Locus/Cell_Type pair is calculated as the simple mean over all expression
#' values measured for that combination.
#'
#' @param ... options to \code{bin_method} and \code{mean_method} (if supplied).
#' Also options for the Efroni binning procedure (\code{l} and \code{u})
#'
#' @return Returns a \code{tibble} containing spec scores and "mean" expression
#' values by Locus. Note that this method does not set negative spec scores to
#' 0. Loci that have an unexpected expression distribution using the Efroni
#' method (background bin is greater than the specified u paramter) will have
#' spec_scores set to 0, however the optimum bin size is still computed, and
#' mean expression is computed using that value.
#'
#' @references
#' Ifroni, E., Ip, PL., Nawy, T., Mello, A., Birnbaum, KD. (2015). "Quantification
#' of cell identity from single-cell gene expression profiles". Genome Biology
#' 16(9)
#'
#' Birnbaum, KD. and Kussell, E. (2011). "Measuring cell identity in noisy
#' biological systems". Nucl. Acids. Res. 39(21)
#'
#' @export
#'
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data
#'
#' @example inst/examples/compute_spec_table_examples.R
#'
compute_spec_table <- function(expression_data,
                               bin_method = "Efroni",
                               mean_method = "Efroni",
                               ...) {

  n_datasets <- unique(expression_data$Sample_Name) %>% length()

  # Pull loci from expression_data structure with complete expresison profiles
  all_loci <- expression_data %>%
    dplyr::filter(!is.na(.data$Expression)) %>%
    dplyr::group_by(.data$Locus) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    dplyr::filter(.data$n == n_datasets) %>%
    dplyr::pull(.data$Locus) %>%
    unique()

  expression_data_filtered <- expression_data %>%
    dplyr::filter(.data$Locus %in% all_loci)

  # For each locus, calculate cell_type mean and specifications scores
  specs <- furrr::future_map_dfr(.x = all_loci, .progress = T, ..., .f = function(x, ...) {
    # Extract locus-specific data
    locus_data <- dplyr::filter(expression_data_filtered, .data$Locus == x)

    # Find optimum bin size (yielding highest spec score)
    # Calculate spec score and mean expression across cell type
    binned_expression_data <- bin_expression_data(
      expression_data = locus_data,
      method = bin_method,
      ...)

    spec_scores <- compute_spec(binned_expression_data)

    mean_expr <- compute_mean_expr(
      binned_expression_data, method = mean_method, ...)

    dplyr::full_join(spec_scores, mean_expr, by = "Cell_Type") %>%
      dplyr::mutate(Locus = x)})

  # Return information in a melted data frame
  return(specs)
}


#' Compute bins for expression data
#'
#' @description takes in expression and cell-type data for a single locus, and
#' computes expression bins.
#'
#' @param expression_data data frame (or tibble) containing the following
#' columns for a single locus
#'
#' \itemize{
#'   \item Expression (the normalized expression value, not log-transformed)
#'   \item Cell_Type (the cell type from which the expression data came from)
#' }
#'
#' @param method a character (possible values consist only of "Efroni" for now)
#' OR a function that takes in the `expression_data` data.frame object and any
#' other parameters, and returns expression_data with an extra column, "bin"
#' containing the expression bin for that dataset/cell_type combination
#' @param ... Additional parameters to the `method` function argument, if specified
#'
#' @return a data.frame (or tibble) containing the data within `expression_data`
#' and an additional "bin" column describing the expression bin for that
#' expression/cell_type combination.
#'
#' @references
#' Ifroni, E., Ip, PL., Nawy, T., Mello, A., Birnbaum, KD. 2015. "Quantification of cell identity from single-cell gene expression profiles". Genome Biology 16(9)
#'
#' Birnbaum, KD. and Kussell, E. 2011. "Measuring cell identity in noisy biological systems". Nucl. Acids. Res. 39(21)
#'
#' @importFrom dplyr "%>%"
bin_expression_data <- function(expression_data, method = "Efroni", ...) {
  if(is.null(method)) method = "NULL"
  if(length(method) == 0) method = "NULL"

  if(class(method) == "function") {
    binned_expression_data <- method(expression_data, ...)
    if(!("bin" %in% colnames(binned_expression_data))) {
      stop(paste0("binning method must return a data.frame with a column, \"bin\""))
    }
  } else if(method == "Efroni") {
    binned_expression_data <- efroni_bin(expression_data, ...)
  } else {
    stop(paste0("bin_method, ",
                method,
                "not recognized (must be a function or \"Efroni\""))
  }
  return(binned_expression_data)
}

#' Comptue mean expression for an `expression_data` object
#'
#' @param expression_data data frame (or tibble) containing the following
#' columns for a single locus
#'
#' \itemize{
#'   \item Expression (the normalized expression value, not log-transformed)
#'   \item Cell_Type (the cell type from which the expression data came from)
#'   \item binsize (optional; required if `method = "Efroni"`)
#' }
#'
#' @param method Which method to use when computing the mean. Defaults to
#' "Efroni", but can be another function. If set to something other than
#' "Efroni", a function, or "median", defaults to computing the median
#' expression for each cell type. Must return one or more summarized expression
#' values for each cell type.
#'
#' @param ... Parameters to be supplied to the function defined in `method`.
#'
#' @return Summarized expression values for each cell type.
#'
#' @references
#' Ifroni, E., Ip, PL., Nawy, T., Mello, A., Birnbaum, KD. 2015. "Quantification of cell identity from single-cell gene expression profiles". Genome Biology 16(9)
#'
#' Birnbaum, KD. and Kussell, E. 2011. "Measuring cell identity in noisy biological systems". Nucl. Acids. Res. 39(21)
#'
#' @importFrom dplyr "%>%"
#' @importFrom stats "median"
#' @importFrom rlang .data
compute_mean_expr <- function(expression_data, method = "Efroni", ...) {
  if(is.null(method)) {
    warning("mean_method set to NULL; setting to \"median\"")
    method = "median"
  }
  if(class(method) == "character") {
    if(!(method %in% c("Efroni", "median"))) {
      warning("mean_method not recognized; setting to \"median\"")
      method = "median"
    }
  }
  if(class(method) == "function") {
    mean_expr <- method(expression_data, ...)
    if(!("mean_expr" %in% colnames(mean_expr))) {
      stop("mean method must return a data.frame with a column, \"mean_expr\"")
    }
  } else if(method == "Efroni") {
    mean_expr <- efroni_mean(expression_data)
  } else {
    if(method != "median") {
      warning(paste0("mean_method, ", method, ", not recognized, defaulting to ",
                   "median expression for all Cell_Type/Locus combinations"))
    }
    mean_expr <- expression_data %>%
      dplyr::group_by(.data$Cell_Type) %>%
      dplyr::summarize(mean_expr = median(.data$Expression))
  }
  return(mean_expr)
}

#' Compute Specification Score
#'
#' @description For expression data for a given locus, calculates the
#' specificity score, as outlined in Birnbaum, et al., 2011. Supports any number
#' of bins, which MUST be specified in the `bin` column of the input object.
#'
#'
#' @param binned_expression_data data frame (or tibble) containing the following
#' columns for a single locus
#'
#' \itemize{
#'   \item Expression (the normalized expression value, not log-transformed)
#'   \item Cell_Type (the cell type from which the expression data came from)
#'   \item bin (which expression bin the expression for this locus falls in per dataset)
#' }
#'
#' @return A tibble with cell-type specific spec scores for this locus. Supports
#' negative spec scores (determined by whether the mean expression for that cell
#' type is above/below the mean expression over all data sets.)
#'
#' @export
#'
#' @references
#' Ifroni, E., Ip, PL., Nawy, T., Mello, A., Birnbaum, KD. 2015. "Quantification of cell identity from single-cell gene expression profiles". Genome Biology 16(9)
#'
#' Birnbaum, KD. and Kussell, E. 2011. "Measuring cell identity in noisy biological systems". Nucl. Acids. Res. 39(21)
#'
compute_spec <- function(binned_expression_data) {
  # Extract information from input data
  mean_expr <- mean(binned_expression_data$Expression)
  ct_vec <- binned_expression_data$Cell_Type
  types <- unique(ct_vec)
  n_types <- length(types)

  bin_table <- table(
    binned_expression_data$Cell_Type,
    binned_expression_data$bin)

  p_level_type <- bin_table/rowSums(bin_table)

  p_type_level <- apply(p_level_type, 2, function(n_bin) n_bin/sum(n_bin))

  information <- apply(p_type_level, 2, function(ptl) 1 + sum(ptl*log(ptl, n_types), na.rm = T))

  spec <- apply(p_level_type, 1, function(plt) sum(information * plt))
  sign <- sign(tapply(binned_expression_data$Expression, binned_expression_data$Cell_Type, mean) - mean_expr)

  return(tibble::enframe(spec*sign, name = "Cell_Type", value = "spec"))
}

#' Convert melted data.frame (or tibble) to matrix
#'
#' @description internal helper function for testing
#'
#' @param df the data.frame
#' @param var1 character, Which variable in df has the new column names
#' @param var2 character, Which variable in df has the new values
#' @param indexvar character, which variable in df has the new rownames
#'
#' @return a matrix
#'
melted_to_mat <- function(df, var1, var2, indexvar) {
  new_df <- dplyr::select(df, var1, var2, indexvar)
  mat <- tidyr::spread(new_df, var1, var2) %>%
    tibble::column_to_rownames(indexvar) %>%
    as.matrix()
}
