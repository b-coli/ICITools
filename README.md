ICITools
========

ICI Tools is an R package implementation of the cell type classification method described in Birnbaum, et al. (2011) and Efroni, et al. (2015). The package implements the computation of a Specificity Score table given a matrix of cell-type specific gene expression data. The package also implements the computation of an Index of Cell Identity score (ICI score) given the specificity score table and new gene expression data.

Installation
------------

To install from github:

    devtools::install_github("b-coli/ICITools")

Computing Spec Table
--------------------

To build a specificity score table using ICITools, you must first create an expression data frame, minimally containing 4 columns:

-   Locus (the gene name)
-   Cell\_Type (which cell type for this gene in this record)
-   Sample\_Name (which sample the expression value came from, can support multiple samples per cell type)
-   Expression (the normalized expression value)

``` r
library(ICITools)
expression_data <- test_spec
head(expression_data)
```

    ## # A tibble: 6 x 4
    ##   Locus   Cell_Type Sample_Name Expression
    ##   <chr>   <chr>     <chr>            <dbl>
    ## 1 locus_1 A         Sample_1        0.0162
    ## 2 locus_2 A         Sample_1        0.690 
    ## 3 locus_3 A         Sample_1        0.229 
    ## 4 locus_4 A         Sample_1        0.141 
    ## 5 locus_5 A         Sample_1        0.0183
    ## 6 locus_6 A         Sample_1        0.0939

A spec table can be computed by running the function, `compute_spec_table()` with default parameters:

``` r
compute_spec_table(expression_data = expression_data)
```

    ## # A tibble: 126 x 6
    ##    Cell_Type    spec binsize mean_expr median_expr Locus   
    ##    <chr>       <dbl>   <dbl>     <dbl>       <dbl> <chr>   
    ##  1 A          0.566    0.205     2.32       0.476  locus_1 
    ##  2 B         -0.0812   0.205     0.449      0.0921 locus_1 
    ##  3 C         -0.0812   0.205     0.323      0.0663 locus_1 
    ##  4 D         -0.0812   0.205     0.292      0.0599 locus_1 
    ##  5 E          0.687    0.205     2.49       0.512  locus_1 
    ##  6 F         -0.0812   0.205     0.498      0.102  locus_1 
    ##  7 G         -0.0812   0.205     0.142      0.0291 locus_1 
    ##  8 H         -0.0812   0.205     0.252      0.0517 locus_1 
    ##  9 I         -0.0812   0.205     0.111      0.0228 locus_1 
    ## 10 A          0        0.102     0.539      0.0551 locus_10
    ## # ... with 116 more rows

This will run the spec score computation using the binning and mean expression computation outlined in Efroni, et al. (2015). To supply your own binning method, define this method as a function, and pass the function into the `bin_method` argument. Note, the function definition must include the `...` parameter to allow for other custom functions in the `compute_spec_table` function.

If a custom bin method is defined, you must also define a custom method to compute expression mean (since the "Efroni" method uses the empirically determined bin size for this). Defining a custom mean function is the same as a custom bin function:

``` r
simple_bin <- function(expression_data, n_bins = 10, ...) {
  bins <- cut(expression_data$Expression, n_bins, labels = FALSE)
  new_expression_data <- dplyr::mutate(expression_data, bin = bins)
  return(new_expression_data)
}

simple_mean <- function(expression_data, ...) {
  means_raw <- tapply(expression_data$Expression, expression_data$Cell_Type, mean)
  means <- tibble::enframe(means_raw, "Cell_Type", "mean_expr")
  return(means)
}

spec_table <- compute_spec_table(expression_data = expression_data, 
                   bin_method = simple_bin, 
                   mean_method = simple_mean, 
                   n_bins = 10)
head(spec_table)
```

    ## # A tibble: 6 x 4
    ##   Cell_Type    spec mean_expr Locus  
    ##   <chr>       <dbl>     <dbl> <chr>  
    ## 1 A          0.690     0.461  locus_1
    ## 2 B         -0.122     0.101  locus_1
    ## 3 C         -0.122     0.0845 locus_1
    ## 4 D         -0.122     0.0674 locus_1
    ## 5 E          0.730     0.595  locus_1
    ## 6 F         -0.0828    0.102  locus_1

The resulting spec table contains spec scores for all loci for all cell types provided. This can be used to compute ICI scores for new data.

Computing ICI scores
--------------------

To compute ICI scores for a new dataset, this dataset must be in a format where each expression profile is an individual column in a data frame, with one column named, "Locus".

``` r
expression_data <- test_ici
head(expression_data)
```

    ## # A tibble: 6 x 16
    ##   Locus Cell_1 Cell_10 Cell_11 Cell_12 Cell_13 Cell_14 Cell_15 Cell_2
    ##   <chr>  <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>  <dbl>
    ## 1 locu~ 5.06    2.28    0.0268   1.29   0.575   0.0713  0.158  0.958 
    ## 2 locu~ 0.155   0.0366  1.23     0.120  1.16    0.0123  0.898  0.112 
    ## 3 locu~ 0.0400  4.29    0.376    0.303  0.0825  0.411   1.13   0.365 
    ## 4 locu~ 0.385   1.86    1.31     0.766  0.322   0.302   0.0522 0.0712
    ## 5 locu~ 0.747   2.64    1.08     0.161  0.352   1.26    1.43   0.183 
    ## 6 locu~ 0.709   2.65    0.413    1.26   0.142   0.115   0.686  0.0335
    ## # ... with 7 more variables: Cell_3 <dbl>, Cell_4 <dbl>, Cell_5 <dbl>,
    ## #   Cell_6 <dbl>, Cell_7 <dbl>, Cell_8 <dbl>, Cell_9 <dbl>

Then, use the function, `compute_ici_scores`. You may specify whether you want significance scores calculated, the amount of information for each cell type (to standardize the amount of specificity conferred by a variable number of markers for each cell type), the number of permutations (for significance computation), and the minimum spec score to be considered a "marker" for each gene.

``` r
ici_scores <- compute_ici_scores(expression_data = expression_data, 
                   spec_table = spec_table, 
                   sig = TRUE, 
                   n_iterations = 1000, 
                   information_level = 20, 
                   min_spec_score = 0.15)
head(ici_scores)
```

    ## # A tibble: 6 x 6
    ##   Cell_Type ici_score p_val ici_score_norm Cell   p_adj
    ##   <chr>         <dbl> <dbl>          <dbl> <chr>  <dbl>
    ## 1 F             0.286 0.988         0.0367 Cell_1 1    
    ## 2 B             1.19  0.626         0.153  Cell_1 0.985
    ## 3 C             0.186 0.968         0.0239 Cell_1 1    
    ## 4 A             3.74  0             0.480  Cell_1 0    
    ## 5 E             1.31  0.566         0.168  Cell_1 0.985
    ## 6 D             0.377 0.798         0.0484 Cell_1 0.985

Parallelization
---------------

This package utilizes the `future` implementation of parallelization for most computationally intensive functions. To enable parallelization, simply call the `future::plan()` function, specifying the number of processors to use. This can consume a significant amount of memory, especially if evaluating a large number of datasets.

``` r
## Time without parallel processing
future::plan(strategy = "sequential")
time <- Sys.time()
ici_scores <- compute_ici_scores(expression_data, spec_table = spec_table, sig = TRUE, n_iterations = 5000)
Sys.time() - time

## Time using 2 processors:
future::plan(strategy = "multiprocess", workers = 2)
time <- Sys.time()
ici_scores <- compute_ici_scores(expression_data, spec_table = spec_table, sig = TRUE, n_iterations = 5000)
Sys.time() - time
```

Optimization
------------

The optimization functions to find an information score that has high signal and low variance, originally outlined in Efroni, et al. (2016), are also available in this package. Here, you can generate ICI scores across a broad range of information scores, specifying a minimum spec score. Additionally, you have the option of sub-sampling within your `expression_data` dataset to speed up optimization (default 100 samples)

``` r
library(dplyr)
library(ggplot2)
info <- gather_information_level_data(expression_data, 
                                      spec_table = spec_table, 
                                      information_range = seq(0, 5, 0.05), 
                                      min_spec_score = 0.15, 
                                      n_samples = 10)

ici_signal <- extract_signal_from_info_data(info)
ici_var <- extract_var_from_info_data(info)

info_summary <- dplyr::inner_join(ici_var, ici_signal)

head(info_summary)
```

ICI Tools Copyright (c) 2019, The Regents of the University of 
California, through Lawrence Berkeley National Laboratory 
(subject to receipt of any required approvals from the U.S.
Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department 
of Energy and the U.S. Government consequently retains certain rights.  As 
such, the U.S. Government has been granted for itself and others acting on 
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.

