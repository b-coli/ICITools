# compute specificity score on test data using Efroni method
# l and u are Efroni-specific hyperparameters that control the expected shape
# of the binned expression data; l = the number of discrete bins for each
# expression profile, and u is the maximum bin that can be considered
# "background". The expression values are then classified as expressed or not
# expressed and spec is calculated on 2 bins.
spec_scores <-
  compute_spec_table(expression_data = test_spec,
                     bin_method = "Efroni", l = 10, u = 3)

# Compute specificity score on test data adapted from Birnbaum, et al.
# (2011), using a custom binning and mean computation method. Make sure to
# supply the "..." argument to both functions unless they use exactly the
# arguments (e.g. if in the below example, ... would not be necessary in
# custom_mean if that function also used the n_bins argument.)

# Custom binning procedure
custom_bin <- function(df, n_bins, ...) {
  bins <- cut(df$Expression, n_bins, labels = FALSE)
  df$bin = bins
  return(df)
}

custom_mean <- function(df, ...) {
  means_raw <- tapply(df$Expression, df$Cell_Type, mean)
  means <- tibble::enframe(means_raw, "Cell_Type", "mean_expr")
  return(means)
}

spec_scores <-
  compute_spec_table(expression_data = test_birnbaum,
                      bin_method = custom_bin, mean_method = custom_mean,
                      n_bins = 3)
