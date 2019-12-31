## Compute ICI scores for a set of test data:

# Train a spec table given initial expression data
spec_table <- compute_spec_table(expression_data = test_spec,
                                 bin_method = "Efroni",
                                 mean_method = "median",
                                 l = 5, u = 3)

# For new cell profiles, compute ICI scores
ici_scores <- compute_ici_scores(expression_data = test_ici,
                                 spec_table = spec_table,
                                 n_iterations = 1000,
                                 sig = TRUE,
                                 information_level = 10,
                                 min_spec_score = 0)
