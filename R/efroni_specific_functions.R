#' Compute mean expression for each cell type (Efroni method)
#'
#' @description Uses the Efroni method to compute mean expression for each cell
#' type present in `expression_data`. Here, it computes the mean expression,
#' divided by the "binsize" used to classify expression.
#'
#' @param expression_data data frame (or tibble) minimally containing the
#' following columns for a single locus
#'
#' \itemize{
#'   \item Expression (the normalized expression value, not log-transformed)
#'   \item Cell_Type (the cell type from which the expression data came from)
#'   \item binsize (the optimized binsize used to compute expression bins)
#' }
#'
#' @return a data frame (or tibble) containing the summarized mean and
#' median expression values for each cell type.
#'
#' @references
#' Ifroni, E., Ip, PL., Nawy, T., Mello, A., Birnbaum, KD. 2015. "Quantification of cell identity from single-cell gene expression profiles". Genome Biology 16(9)
#'
#' Birnbaum, KD. and Kussell, E. 2011. "Measuring cell identity in noisy biological systems". Nucl. Acids. Res. 39(21)
#'
#' @importFrom dplyr "%>%"
#' @importFrom stats "median"
#' @importFrom rlang .data
efroni_mean <- function(expression_data) {
  expression_data %>%
    dplyr::group_by(.data$Cell_Type) %>%
    dplyr::summarize(
      binsize = unique(.data$binsize),
      mean_expr = dplyr::if_else(.data$binsize > 0, median(.data$Expression/.data$binsize), 0),
      median_expr = median(.data$Expression))
}

#' Optimize the selection of bin size for defining two expression classes
#' in spec score computation.
#'
#' @param expression_data  data frame (or tibble) containing the following
#' columns for a single locus
#'
#' \itemize{
#'   \item Expression (the normalized expression value, not log-transformed)
#'   \item Cell_Type (the cell type from which the expression data came from)
#' }
#'
#' @param l integer, the number of cuts to perform when binning expression values (see Efroni, et al., 2015)
#' @param u integer, the maximum bin to consider as a background bin (0 < u < l)
#' @param use_original logical, whether to use the original method of computing
#' background bins (ignores minimum values). FALSE by default.
#'
#' @references
#' Ifroni, E., Ip, PL., Nawy, T., Mello, A., Birnbaum, KD. 2015. "Quantification
#' of cell identity from single-cell gene expression profiles". Genome Biology
#' 16(9)
#'
#' Birnbaum, KD. and Kussell, E. 2011. "Measuring cell identity in noisy
#' biological systems". Nucl. Acids. Res. 39(21)
#'
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data
efroni_bin <- function(expression_data, l = 10, u = 3, use_original = F) {
  expression <- expression_data$Expression
  n_datasets <- length(expression)
  expr_thresh <- ceiling(n_datasets/l)

  background_bin <- efroni_get_background_bin(expression, expr_thresh, l, use_original)
  binsizes <- efroni_get_binsizes(expression, background_bin, l)

  spec_by_binsize <- purrr::map_dfr(.x = binsizes, .f = function(binsize) {
    binned_expression_data <- data.frame(
      Cell_Type = expression_data$Cell_Type,
      Expression = expression_data$Expression,
      bin = cut(expression_data$Expression, breaks = c(-1, binsize, Inf), labels = F))
    spec_raw <- compute_spec(binned_expression_data)
    return(cbind(spec_raw, binsize))
  })

  optimum_binsize <- dplyr::group_by(spec_by_binsize, .data$binsize) %>%
    dplyr::summarize(max_spec = max(.data$spec, na.rm = T)) %>%
    dplyr::top_n(1, .data$max_spec) %>%
    dplyr::pull(.data$binsize) %>% min(na.rm = T)

  if(!is.na(background_bin) & background_bin > u) {
    expression_breaks = c(min(expression) - 1, max(expression) + 1)
  } else {
    expression_breaks = c(min(expression)-1, optimum_binsize, Inf)
  }

  binned_expression_data <- expression_data %>%
    dplyr::mutate(binsize = optimum_binsize) %>%
    dplyr::mutate(bin = cut(.data$Expression,
                            breaks = expression_breaks,
                            labels = F))
  return(binned_expression_data)
}

#' Find the background bin in an expression vector
#'
#' @description Implements modified Efroni method for determining which bin
#' corresponds to background expression. This method is used internally within
#' the `bin_expression_data`` function.
#'
#' @param expression a numeric vector of expression values
#' @param expression_thresh the threhold for determining whether a bin is in the background
#' @param l integer, number of expression bins
#' @param use_original logical, whether to use the original method of computing
#' background bins (ignores minimum values). FALSE by default.
#'
#' @return an integer corresponding to the "background bin". If this bin is
#' undetermined (no clear grouping of data), returns the maximum bin number (l)
#'
#' @references
#' Ifroni, E., Ip, PL., Nawy, T., Mello, A., Birnbaum, KD. 2015. "Quantification
#' of cell identity from single-cell gene expression profiles". Genome Biology
#' 16(9)
#'
#' Birnbaum, KD. and Kussell, E. 2011. "Measuring cell identity in noisy
#' biological systems". Nucl. Acids. Res. 39(21)
#'
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data
efroni_get_background_bin <- function(expression, expression_thresh, l, use_original = F) {
  if(l <= 1) return(NA)

  if(use_original) {
    expression_range <- max(expression)-min(expression)
    expression_breaks <- seq(min(expression),max(expression),expression_range/l)
    expression_bin_table <- cut(expression, breaks = expression_breaks, labels = F) %>%
      table() %>% tibble::enframe("bin", "n") %>%
      dplyr::mutate(bin = as.numeric(.data$bin))
  } else {
    expression_bin_table <- cut(expression, breaks = l, labels = F) %>%
      table() %>% tibble::enframe("bin", "n") %>%
      dplyr::mutate(bin = as.numeric(.data$bin))
  }

  background_bin <- expression_bin_table %>%
    dplyr::filter(.data$n >= expression_thresh) %>%
    dplyr::top_n(1, .data$bin) %>%
    dplyr::pull(.data$bin) %>%
    as.numeric()

  if(length(background_bin) > 0) {
    return(background_bin)
  } else {
    return(l)
  }
}

#' Get binsize range for optimization
#'
#' @description If the number of cuts to test are specified, compute the range
#' of possible bin sizes as close to the transition to "background" expression
#' as possible. If the number of cuts is not specified, attempts to optimize the binsize across the whole range of gene expression values (50 evenly spaced binsizes).
#'
#' @param expression numeric, expression vector
#' @param background_bin pre-computed background bin. If NA, whole expression
#' range is considered.
#' @param l number of cuts (bins). If < 1, whole expression range is considered.
#'
#' @return a vector of 50 bin sizes to iterate over.
#'
#' @references
#' Ifroni, E., Ip, PL., Nawy, T., Mello, A., Birnbaum, KD. 2015. "Quantification
#' of cell identity from single-cell gene expression profiles". Genome Biology
#' 16(9)
#'
#' Birnbaum, KD. and Kussell, E. 2011. "Measuring cell identity in noisy
#' biological systems". Nucl. Acids. Res. 39(21)
#'
efroni_get_binsizes <- function(expression, background_bin, l) {

  expression_range <- max(expression)-min(expression)

  if(l > 0 & !is.na(background_bin)) {
    binsizes_min <- min(expression) + ((expression_range/l)*(background_bin - 1.5))
    binsizes_max <- min(expression) + ((expression_range/l)*(background_bin + 0.5))
  } else {
    binsizes_min <- min(expression) + (expression_range)*0.1
    binsizes_max <- max(expression) - (expression_range)*0.1
  }

  binsizes_range <- binsizes_max - binsizes_min
  binsizes <- seq(binsizes_min, binsizes_max, binsizes_range/50)

  return(binsizes)
}
