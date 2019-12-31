test_that("ici computation works", {
  spec_table <- compute_spec_table(test_spec, l = 5, u = 2)
  ici_scores <- compute_ici_scores(expression_data = test_ici,
                                   spec_table = spec_table,
                                   information_level = 10,
                                   min_spec_score = 0,
                                   sig = TRUE)
  top_ici_scores <- ici_scores %>%
    dplyr::group_by(Cell) %>%
    dplyr::top_n(1, ici_score_norm) %>%
    dplyr::select(Cell, Cell_Type) %>%
    tidyr::spread(Cell, Cell_Type)

  expect_equal(top_ici_scores$Cell_1, "A")
  expect_equal(top_ici_scores$Cell_2, "B")
  expect_equal(top_ici_scores$Cell_3, "C")
  expect_equal(top_ici_scores$Cell_4, "D")
  expect_equal(top_ici_scores$Cell_5, "E")
  expect_equal(top_ici_scores$Cell_6, "F")
})


test_that("compute_ici returns p-values for expression data missing some loci", {
  test_ici_missing <- dplyr::filter(test_ici, Locus != "locus_4")
  spec_table <- compute_spec_table(test_spec, l = 5, u = 3)
  test_ici_missing_scores <- compute_ici_scores(expression_data = test_ici_missing,
                                                spec_table = spec_table,
                                                sig = T, information_level = 4,
                                                min_spec_score = 0)
  n_missing <- sum(is.na(test_ici_missing_scores$p_val))
  expect_equal(n_missing, 0)
})
