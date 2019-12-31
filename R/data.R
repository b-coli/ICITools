#' Toy dataset adapted from Birnbaum, et al., 2011
#'
#' @format A dataset of 1 locus, measured across 6 hypothetical cell types (2-4 replicates per cell type).
#' Conforms to the expected expression_data format, with 3 columns:
#'
#' \describe{
#'   \item{Locus}{Locus identifer (chr)}
#'   \item{Expression}{Expression of this locus in a given dataset (dbl)}
#'   \item{Cell_Type}{cell type of the data set (chr)}
#'   \item{Sample_Name}{name of the "sample" this locus/Cell_Type measurement originated from}
#' }
#'
#' Can be used to test the compute_spec_table function. To reproduce results
#' reported in Birnbaum, et al., (2011), supply the function:
#' `bin_func <- function(df) dplyr::mutate(df, bin = cut(Expression, 3))`
#' as the bin_method, and null string, "", as the mean_method to
#' compute_spec_table.
"test_birnbaum"

#' Random toy dataset simulating marker loci
#'
#' @format A dataset of 10 loci, measured across 5 hypothetical cell types (5
#' replicates per cell type). Conforms to the expected expression_data format,
#' with 4 columns:
#'
#' \describe{
#'   \item{Locus}{Locus identifer (chr)}
#'   \item{Expression}{Expression of this locus in a given dataset (dbl)}
#'   \item{Cell_Type}{cell type of the data set (chr)}
#'   \item{Sample_Name}{name of the "sample" this locus/Cell_Type measurement originated from}
#' }
#'
#' Can be used to test the compute_spec_table function.
"test_spec"

#' Random toy dataset simulating cell profiles
#'
#' @format A dataset of 15 cells, with the same 14 loci as described in
#' \code{test_spec}. Conforms to the expected expression_data format, with 3
#' columns:
#'
#' \describe{
#'   \item{Locus}{Locus identifer (chr)}
#'   \item{Expression}{Expression of this locus in a given dataset (dbl)}
#'   \item{Cell}{cell type of the data set (chr)}
#' }
#'
#' Can be used to test the compute_ici_scores function.
"test_ici"

