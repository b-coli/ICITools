library(dplyr)

alt_bin_method <- function(df, n_bins, ...) {
  bins <- cut(df$Expression, n_bins, labels = F)
  df$bin = bins
  return(df)
}

test_that("compute_spec_table runs cleanly on test data", {
  result <- compute_spec_table(test_spec)
  expect_true(sum(abs(result$spec)) > 3)
  expect_identical(dim(result), as.integer(c(126,6)))
})

test_that("Birnbaum example with custom bin method is accurate", {
  result <- compute_spec_table(test_birnbaum,
                               bin_method = alt_bin_method,
                               mean_method = "median",
                               n_bins = 3)
  expected_spec <- c(0.12, 0.29, 0.72, 0.12, 0.40, 0.12)
  expect_true(cor(abs(result$spec), expected_spec) > 0.99)
})

test_that("future works on data sets", {
  future::plan(strategy = "sequential")
  time <- Sys.time()
  for(i in 1:10) result <- compute_spec_table(test_spec)
  t1 <- Sys.time() - time

  future::plan(strategy = "multiprocess")
  time <- Sys.time()
  for(i in 1:10) result <- compute_spec_table(test_spec)
  t2 <- Sys.time() - time

  result <- sum(abs(result$spec))
  if(future::availableCores() > 1) {
    expect_true(result > 0)
    expect_true(t2 < t1)
  } else {
    expect_true(result > 0)
  }
})

test_that("compute_spec_table handles other inputs to bin_method and mean_method", {
  expect_error(compute_spec_table(test_spec, bin_method = "my_method"))
  expect_error(compute_spec_table(test_spec, bin_method = NULL))
  expect_error(compute_spec_table(test_spec, bin_method = function(x) {}))

  expect_warning(compute_spec_table(test_spec, mean_method = "my_method"))
  expect_warning(compute_spec_table(test_spec, mean_method = NULL))
  expect_error(compute_spec_table(test_spec, mean_method = function(x) {}))

  expect_error(compute_spec_table(test_spec, bin_method = "my_method", mean_method = "my_method"))
  expect_error(compute_spec_table(test_spec, bin_method = NULL, mean_method = NULL))
  expect_error(compute_spec_table(test_spec, bin_method = function(x) {}, mean_method = function(x) {}))
})
