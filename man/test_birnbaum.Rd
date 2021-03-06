% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{test_birnbaum}
\alias{test_birnbaum}
\title{Toy dataset adapted from Birnbaum, et al., 2011}
\format{A dataset of 1 locus, measured across 6 hypothetical cell types (2-4 replicates per cell type).
Conforms to the expected expression_data format, with 3 columns:

\describe{
  \item{Locus}{Locus identifer (chr)}
  \item{Expression}{Expression of this locus in a given dataset (dbl)}
  \item{Cell_Type}{cell type of the data set (chr)}
  \item{Sample_Name}{name of the "sample" this locus/Cell_Type measurement originated from}
}

Can be used to test the compute_spec_table function. To reproduce results
reported in Birnbaum, et al., (2011), supply the function:
`bin_func <- function(df) dplyr::mutate(df, bin = cut(Expression, 3))`
as the bin_method, and null string, "", as the mean_method to
compute_spec_table.}
\usage{
test_birnbaum
}
\description{
Toy dataset adapted from Birnbaum, et al., 2011
}
\references{
Birnbaum, KD. and Kussell, E. 2011. "Measuring cell identity in noisy
biological systems". Nucl. Acids. Res. 39(21)
}
\keyword{datasets}
