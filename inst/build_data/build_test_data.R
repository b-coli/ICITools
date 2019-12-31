test_spec <- tibble::tibble(
  Locus = rep(paste0("locus_", 1:14), 35),
  Cell_Type = c(rep(LETTERS[1:5], 5, each = 14),
                rep("F", 1, each = 14),
                rep("G", 2, each = 14),
                rep("H", 3, each = 14),
                rep("I", 4, each = 14)),
  Sample_Name = paste0("Sample_", rep(1:35, each = 14))) %>%
  dplyr::mutate(
    Expression = dplyr::case_when(
      Cell_Type == "A" & Locus %in% c("locus_1", "locus_2", "locus_3", "locus_4") ~ 1,
      Cell_Type == "B" & Locus %in% c("locus_2", "locus_5", "locus_6") ~ 1,
      Cell_Type == "C" & Locus %in% c("locus_7") ~ 1,
      Cell_Type == "D" & Locus %in% c("locus_8", "locus_9") ~ 1,
      Cell_Type == "E" & Locus %in% c("locus_1", "locus_9") ~ 2,
      Cell_Type == "F" & Locus %in% c("locus_10", "locus_11", "locus_12") ~ 1
    )
  ) %>% dplyr::mutate(Expression = dplyr::if_else(is.na(Expression), 0, Expression))

set.seed(1)
test_spec$Expression = test_spec$Expression + runif(length(test_spec$Expression), 0, 0.35)
test_spec$Expression = test_spec$Expression*runif(length(test_spec$Expression), 0, 0.65)
usethis::use_data(test_spec, internal = F, overwrite = T)

test_birnbaum <- tibble::tibble(
  Locus = "locus_1",
  Cell_Type = c(rep(c("A", "B", "C", "D", "E"), each = 3), "E", rep("F", 2)),
  Expression = c(0, 0.2, 0.5,
                 2, 0.4, 0.2,
                 4, 5, 4.5,
                 0.1, 0, 0.2,
                 5, 2.3, 0.3,
                 0.1, 0.1, 0),
  Sample_Name = paste0("Sample_", 1:18)
)
usethis::use_data(test_birnbaum, internal = F, overwrite = T)

test_ici <- tibble::tibble(
  Locus = rep(paste0("locus_", 1:14), 15),
  Cell = rep(paste0("Cell_", 1:15), each = 14)) %>%
  dplyr::mutate(Expression = dplyr::case_when(
    Cell == "Cell_1" & Locus %in% paste0("locus_", 1:4) ~ 50,
    Cell == "Cell_2" & Locus %in% paste0("locus_", c(2,5,6)) ~ 50,
    Cell == "Cell_3" & Locus == "locus_7" ~ 75,
    Cell == "Cell_4" & Locus %in% paste0("locus_", c(8,9)) ~ 50,
    Cell == "Cell_5" & Locus %in% paste0("locus_", c(1,9)) ~ 200,
    Cell == "Cell_6" & Locus %in% paste0("locus_", c(10,11,12)) ~ 25,
    Cell == "Cell_7" & Locus == "locus_14" ~ 35,
    Cell == "Cell_8" & Locus == "locus_8" ~ 91,
    Cell == "Cell_9" & Locus == "locus_1" ~ 9,
    Cell == "Cell_10" & Locus != "locus_10" ~ 16
  )) %>% dplyr::mutate(Expression = dplyr::if_else(is.na(Expression), 0, Expression))
set.seed(15)
test_ici$Expression = test_ici$Expression + runif(length(test_ici$Expression), 0, 8)
test_ici$Expression = test_ici$Expression*runif(length(test_ici$Expression), 0, 0.25)
test_ici <- tidyr::spread(test_ici, Cell, Expression)
usethis::use_data(test_ici, internal = F, overwrite = T)

download.file(url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4354993/bin/13059_2015_580_MOESM5_ESM.zip",
              destfile = "inst/efroni_implementation.zip")
unzip("inst/efroni_implementation.zip", exdir = "inst")
source("inst/identity.r")
source("inst/spec.R")
