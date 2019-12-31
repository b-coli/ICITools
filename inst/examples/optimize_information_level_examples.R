# Try to find the optimum information level
spec_table <- compute_spec_table(test_spec, l = 5, u = 3)
information_level <- optimize_information_level(expression_data = test_ici,
                                                spec_table = spec_table,
                                                information_range = seq(0.5,5,0.5),
                                                min_spec_score = 0,
                                                n_samples = 15,
                                                var_tol = 0.2)

# Gather information for all information levels
info_data <- gather_information_level_data(test_ici,
                                           spec_table,
                                           information_range = seq(0.5,5,0.5),
                                           min_spec_score = 0,
                                           n_samples = 15)
