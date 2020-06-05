#' Find Markers for ICI score computation
#'
#' @param spec_table a data.frame containing at least 4 columns:
#' \itemize{
#'   \item Locus
#'   \item Cell_Type
#'   \item spec (the Specificity score calculated using compute_spec_table)
#'   \item mean_expr (the mean expression level for that locus/cell type
#'   combination used to generated spec scores)
#' }
#'
#' @param information_level numeric, the maximum cumulative information level
#' used to compute the ICI score
#'
#' @param transform_spec logical, whether to convert spec scores below
#' min_spec_score to 0 (default: TRUE)
#'
#' @param min_spec_score numeric, minimum useful specification score
#'
#' @return Transforms the input spec_table with two additional columns:
#' \itemize{
#'   \item cumulative_sum, the cumulative sum of spec scores descending from the
#'   highest to lowest (with mean_expr used to break ties)
#'   \item is_marker whether the cumulative sum for that locus is below the
#'   information threshold specified by information_level
#' }
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
find_ici_markers <- function(spec_table, information_level, min_spec_score = 0.15, transform_spec = TRUE) {
  if(transform_spec) {
    spec_table <- dplyr::mutate(
      spec_table,
      spec = dplyr::if_else(.data$spec > min_spec_score, .data$spec, 0))
  }
  spec_table %>%
    dplyr::group_by(.data$Cell_Type) %>% # Loop over all cell types
    dplyr::arrange(dplyr::desc(.data$spec), dplyr::desc(.data$mean_expr)) %>% # Order spec scores by decreasing spec and mean expression
    dplyr::mutate(cumulative_sum = cumsum(.data$spec)) %>% # Compute cumulative spec score
    dplyr::mutate(is_marker = .data$cumulative_sum < information_level) %>% # Filter loci with spec scores lower than information level
    dplyr::mutate(is_marker = dplyr::if_else(.data$spec > min_spec_score, .data$is_marker, FALSE)) %>%
    dplyr::ungroup() # Flatten data frame
}

#' Compute ICI score for a given cell type (Internal)
#'
#' @param marker_expression numeric vector, gene expression profile
#' @param marker_spec numeric vector, specification scores for corresponding genes in marker_expression
#'
#' @return numeric, index of cell identity score for this expression/spec score combination
#'
compute_ici_score_for_cell_type <- function(marker_expression, marker_spec) {
  score.1 <- mean(marker_expression * marker_spec)
  score.2 <- mean(marker_expression > 0)
  return(score.1 * score.2)
}

#' Compute p-values for a single expression/cell type profile (Internal)
#'
#' @param cell_type_data a data.frame with the following columns:
#' \itemize{
#'   \item Locus, the gene name
#'   \item is_marker, logical, whether each gene is a marker
#'   \item Expression, numeric, the expression of each locus for this profile
#'   \item spec, numeric, the specificity score for this each locus for the cell
#'   type in question.
#'  }
#' @param ici_score numeric, the ici score to compare the random iterations against.
#' If null (default), will compute denovo
#'
#' @param n_iterations, integer, number of random iterations
#'
#' @return numeric, raw p-value indicating the number of randomized ici scores
#' that are greater than or equal to the "real" ici score.
#' 
#' @references
#' Ifroni, E., Ip, PL., Nawy, T., Mello, A., Birnbaum, KD. 2015. "Quantification
#' of cell identity from single-cell gene expression profiles". Genome Biology
#' 16(9)
#'
#' Birnbaum, KD. and Kussell, E. 2011. "Measuring cell identity in noisy
#' biological systems". Nucl. Acids. Res. 39(21)
#'
#'
#' @importFrom dplyr "%>%"
compute_ici_pval <- function(cell_type_data, ici_score = NULL, n_iterations = 1000) {
  if(is.null(ici_score)) {
    ici_score = compute_ici_score_for_cell_type(
      marker_expression = cell_type_data$Expression,
      marker_spec = cell_type_data$spec
    )
  }

  n_markers <- sum(cell_type_data$is_marker)
  n_genes <- length(cell_type_data$Locus)
  marker_index <- which(cell_type_data$is_marker)

  set.seed(1)
  scores <- purrr::map(1:n_iterations, .f = function(x) {
    random_index <- sample(1:n_genes, n_markers)
    random_marker_expression <- cell_type_data$Expression[random_index]
    random_marker_spec <- cell_type_data$spec[marker_index]
    spec <- compute_ici_score_for_cell_type(random_marker_expression, random_marker_spec)
    return(spec)}) %>% purrr::reduce(c)
  p <- sum(scores >= ici_score)/n_iterations
}


#' Compute ici score for all cell types for a given cell's expression profile
#'
#' @param expression_and_marker_data data.frame with the following columns for
#' each cell type (melted):
#' \itemize{
#'   \item Locus, the gene name
#'   \item Cell_Type, the cell type for each Locus/Expression/spec/is_marker combination
#'   \item is_marker, logical, whether each gene is a marker
#'   \item Expression, numeric, the expression of each locus for this profile
#'   \item spec, numeric, the specificity score for this each locus for the cell
#'   type in question.
#'  }
#' @param sig logical, whether to compute p-values for ici scores
#' @param n_iterations integer, how many random permutations for p-value calculations
#'
#' @return a data.frame with three columns:
#' \itemize{
#'   \item Cell_Type, the cell type for each Locus/Expression/spec/is_marker combination
#'   type in question.
#'   \item ici_score the ici score for that cell/cell-type
#'   \item p-val the p-value associated with the ici-score, NA if sig = FALSE.
#'  }
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
compute_ici_for_cell <- function(expression_and_marker_data, sig = FALSE, n_iterations = 1000) {
  cell_types <- unique(expression_and_marker_data$Cell_Type)
  ici_score <- purrr::map_dfr(cell_types, .f = function(cell_type) {
    cell_type_data <- dplyr::filter(expression_and_marker_data, .data$Cell_Type == cell_type)
    marker_data <- dplyr::filter(cell_type_data, .data$is_marker)
    ici_score <- compute_ici_score_for_cell_type(marker_data$Expression, marker_data$spec)
    p_val <- NA
    if(sig) p_val <- compute_ici_pval(cell_type_data, ici_score, n_iterations)
    return(tibble::tibble(Cell_Type = cell_type, ici_score, p_val))
  })

  ici_score_processed <- ici_score %>%
    dplyr::mutate(ici_score_norm = .data$ici_score/sum(.data$ici_score)) %>%
    dplyr::mutate(ici_score = dplyr::if_else(is.na(.data$ici_score), 0, .data$ici_score)) %>%
    dplyr::mutate(ici_score_norm = dplyr::if_else(is.na(.data$ici_score_norm), 0, .data$ici_score_norm))
}


#' Compute ICI scores for a set of cells
#'
#' @param expression_data data frame with at least one column, "Locus", and
#' additional columns containing expression values for each cell.
#'
#' @param spec_table a data.frame containing at least 4 columns:
#' \itemize{
#'   \item Locus
#'   \item Cell_Type
#'   \item spec (the Specificity score calculated using compute_spec_table)
#'   \item mean_expr (the mean expression level for that locus/cell type
#'   combination used to generated spec scores)
#' }
#'
#' @param sig logical, whether to compute significance scores for ici values
#' @param n_iterations integer, number of iterations for empirical p-value computation
#' @param information_level numeric, cumulative information score when determining marker loci
#' @param min_spec_score minimum useful spec score for determining marker loci
#'
#' @return a data frame containing ici values and p-values for all cell/cell-type combinations.
#' @example inst/examples/compute_ici_scores_examples.R
#' @export
#'
#' @references
#' Ifroni, E., Ip, PL., Nawy, T., Mello, A., Birnbaum, KD. 2015. "Quantification
#' of cell identity from single-cell gene expression profiles". Genome Biology
#' 16(9)
#'
#' Birnbaum, KD. and Kussell, E. 2011. "Measuring cell identity in noisy
#' biological systems". Nucl. Acids. Res. 39(21)
#'
#' @importFrom stats "p.adjust"
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data
compute_ici_scores <- function(
  expression_data,
  spec_table,
  sig = FALSE,
  n_iterations = 1000,
  information_level = 20,
  min_spec_score = 0.15)
{
  all_cells <- expression_data %>% dplyr::select(-.data$Locus) %>% colnames()
  marker_data <- find_ici_markers(
    spec_table = spec_table,
    information_level = information_level,
    min_spec_score = min_spec_score,
    transform_spec = TRUE)

  expression_data_df <- expression_data %>%
    tibble::column_to_rownames("Locus")

  expression_data_list <- purrr::imap(expression_data_df, function(x, y) {
    df <- as.data.frame(x)
    colnames(df) <- y
    row.names(df) <- row.names(expression_data_df)
    Matrix::Matrix(as.matrix(df), sparse = T)})

  required_pkgs <- c("Matrix")
  f_opts <- furrr::future_options(packages = required_pkgs)

  ici_scores <- furrr::future_map_dfr(
    .x = expression_data_list,
    .options = f_opts,
    .progress = TRUE,
    .f = function(cell_data) {
      cell_name <- colnames(cell_data)
      cell_expression_data <- as.matrix(cell_data) %>%
        tibble::as_tibble(rownames = "Locus") %>%
        `colnames<-`(c("Locus", "Expression"))
      expression_and_marker_data <- dplyr::inner_join(marker_data, cell_expression_data, by = "Locus")
      scores <- compute_ici_for_cell(expression_and_marker_data, sig = sig, n_iterations = n_iterations)
      dplyr::mutate(scores, Cell = cell_name)
    })
  ici_scores <- dplyr::mutate(ici_scores, p_adj = p.adjust(.data$p_val, method = "BH"))
}

#' Attempt to optimize information level
#'
#' @description The information level information to the ICI computation
#' algorithm is curcial for the selection of markers. This function attempts to
#' find an ideal information level by computing the "variability" of ICI scores
#' between sequential information scores, controlling this level at a user-
#' defined rate (var_tol), and then finding the information level that maximizes
#' the ICI signal under this variability. This function is experimental.
#'
#' @param expression_data data frame with at least one column, "Locus", and
#' additional columns containing expression values for each cell.
#'
#' @param spec_table a data.frame containing at least 4 columns:
#' \itemize{
#'   \item Locus
#'   \item Cell_Type
#'   \item spec (the Specificity score calculated using compute_spec_table)
#'   \item mean_expr (the mean expression level for that locus/cell type
#'   combination used to generated spec scores)
#' }
#'
#' @param information_range numeric vector, range of information scores to loop over
#' @param min_spec_score minimum useful information score (default, 0.15)
#' @param n_samples number of sub-samples to take (useful for limiting computational resources required; default 100)
#' @param var_tol numeric (between 0 and 1), the percentile of variation. After
#' computing spec scores for the sampled cells for all cell types at all
#' information levels supplied in information_range, this function then computes
#' the variation (euclidean distance of normalized spec scores over all cell
#' types for a given cell between information levels that differ by 1).
#' All *unique* variation values are then determined among all the cells, and
#' sorted. The value in var_tol is then multiplied by the length of all unique
#' values to determine the variation threshold. Tthe information levels that
#' have maximum variation scores below this threshold are then used to maximize
#' the ICI signal. The maximum information level with the maximum ICI score and
#' variation below this cutoff is returned.
#'
#' @return numeric, information level that should maximize ICI signal
#'
#' @export
#'
#' @example inst/examples/optimize_information_level_examples.R
#'
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data
optimize_information_level <- function(expression_data,
                                       spec_table,
                                       information_range = seq(0,100,0.5),
                                       min_spec_score = 0.15,
                                       n_samples = 100, var_tol = 0.2) {
  info_data <- gather_information_level_data(
    expression_data = expression_data,
    spec_table = spec_table,
    information_range = information_range,
    min_spec_score = min_spec_score,
    n_samples = n_samples)

  var_data <- extract_var_from_info_data(info_data)
  var_vector <- unique(var_data$variability) %>% sort()
  max_var <- max(var_vector[cut(var_vector, 100, labels = F)/100 < var_tol])
  if(is.finite(max_var)) {
    signal_data <- extract_signal_from_info_data(info_data)
    allowed_information_levels <- var_data %>%
      dplyr::group_by(.data$information_level) %>%
      dplyr::summarize(info_max_var = max(.data$variability)) %>%
      dplyr::filter(.data$info_max_var < max_var) %>%
      dplyr::pull(.data$information_level) %>%
      unique()
    optimum_info <- signal_data %>%
      dplyr::filter(.data$information_level %in% allowed_information_levels) %>%
      dplyr::top_n(1, .data$signal) %>%
      dplyr::top_n(1, dplyr::desc(.data$information_level)) %>%
      dplyr::pull(.data$information_level)
    return(optimum_info)
  } else {
    return(NA)
  }
}

#' Extract mean, maximimum ICI signal per cell per information level (internal)
#'
#' @param info_data information data structure returned by gather_information_level_data
#'
#' @importFrom rlang .data
#'
#' @return a data frame with information level/ICI signal data
#'
#' @export
extract_signal_from_info_data <- function(info_data) {
  signal_data <- info_data %>%
    dplyr::group_by(.data$Cell, .data$information_level) %>%
    dplyr::summarize(max_ici = max(.data$ici_score_norm)) %>%
    dplyr::group_by(.data$information_level) %>%
    dplyr::summarize(signal = mean(.data$max_ici))
}

#' Extract variation between information levels for all cells (internal)
#'
#' @param info_data information data structure returned by gather_information_level_data
#'
#' @return a data frame with information level/ICI variation information
#' @importFrom rlang .data
#' @importFrom stats dist
#' @export
extract_var_from_info_data <- function(info_data) {
  var_results <- purrr::map_dfr(unique(info_data$Cell), .f = function(x) {
    cell_info_data <-
      dplyr::filter(info_data, .data$Cell == x) %>%
      dplyr::select(.data$information_level,
                    .data$Cell_Type,
                    .data$ici_score_norm)

      cell_info_data_mat <-
        tidyr::spread(cell_info_data,
                      .data$Cell_Type,
                      .data$ici_score_norm) %>%
        tibble::column_to_rownames("information_level")

      cell_info_dist <- dist(cell_info_data_mat) %>% as.matrix() %>% as.data.frame()

      results <- tibble::rownames_to_column(cell_info_dist, "information_level_1") %>%
        tidyr::gather("information_level_2", "distance", -.data$information_level_1) %>%
        dplyr::filter(as.numeric(.data$information_level_2) == (as.numeric(.data$information_level_1) + 1)) %>%
        dplyr::mutate(information_level = as.numeric(.data$information_level_1)) %>%
        dplyr::select(.data$information_level, variability = .data$distance) %>%
        dplyr::mutate(Cell = x)
      return(results)
  })

  return(var_results)
}

#' Find ICI scores over a range of information levels
#'
#' @param expression_data data frame with at least one column, "Locus", and
#' additional columns containing expression values for each cell.
#'
#' @param spec_table a data.frame containing at least 4 columns:
#' \itemize{
#'   \item Locus
#'   \item Cell_Type
#'   \item spec (the Specificity score calculated using compute_spec_table)
#'   \item mean_expr (the mean expression level for that locus/cell type
#'   combination used to generated spec scores)
#' }
#'
#' @param information_range numeric vector, range of information scores to loop over
#' @param min_spec_score minimum useful information score (default, 0.15)
#' @param n_samples number of sub-samples to take (useful for limiting computational resources required; default 100)
#'
#' @return a data frame with ici scores for all sampled cells, all cell types,
#' at all information levels in information_range.
#'
#' @export
#'
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data
gather_information_level_data <- function(
  expression_data,
  spec_table,
  information_range = seq(0,100,0.5),
  min_spec_score = 0.15,
  n_samples = 100)
{
  all_cells <- colnames(expression_data) %>% setdiff("Locus")
  subset_cells <- sample(all_cells, n_samples)

  expression_data_subset <- dplyr::select(expression_data, .data$Locus, sample(all_cells, n_samples))

  info_data <- furrr::future_map_dfr(
    .progress = TRUE,
    .x = information_range,
    .f = function(information_level) {
      ici_scores <- compute_ici_scores(
        expression_data = expression_data_subset,
        min_spec_score = min_spec_score,
        spec_table = spec_table,
        information_level = information_level,
        sig = FALSE)

      dplyr::mutate(ici_scores, information_level = information_level)
    })
}
