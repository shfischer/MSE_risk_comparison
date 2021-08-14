Exploratory MSE for ple.27.7e - 2021
================

## Introduction

This repository contains the scripts for an exploratory Management
Strategy Evaluation for Western English Channel plaice (ple.27.7e).

The operating model is created using the SAM
[`stockassessment`](https://github.com/fishfollower/SAM/) R package and
follows the principles developed during the ICES Workshop on North Sea
stocks management strategy evaluation
([WKNSMSE](https://doi.org/10.17895/ices.pub.5090)).

The simulation is based on the Fisheries Library in R
([FLR](http://www.flr-project.org/)) and the Assessment for All (a4a)
standard MSE framework ([`FLR/mse`](github.com/FLR/mse)) developed
during the Workshop on development of MSE algorithms with R/FLR/a4a
([Jardim et al.,
2017](https://ec.europa.eu/jrc/en/publication/assessment-all-initiativea4a-workshop-development-mse-algorithms-rflra4a)).

The data-limited MSE framework is based on Fischer et
al. ([2021](https://dx.doi.org/10.1093/icesjms/fsab018)), see the
[`shfischer/GA_MSE`](https://github.com/shfischer/GA_MSE) GitHub
repository for the original source code and documentation.

## Repository structure

The root folder contains the following R scripts:

-   `OM*.R`: Scripts for creating the operating model (OM)
    -   `OM_SAM_fit.R` script for configuring and running SAM,
    -   `OM_Eqsim.R` runs EqSim on SAM model fit to create ICES style
        reference points,
    -   `OM_length.R` preparation of length data for the generation of
        the length index for the rfb rule,
    -   `OM.R` generate the operating model(s),
    -   `OM_MSY.R` calculation of real MSY reference points with the mse
        framework,
-   `funs.R` contains functions and methods used for the creation of the
    operating models and the modules of the MSE,
-   `funs_GA.R` functions for running the MSE within the genetic
    algorithm and calculate summary statistics and fitness values,
-   `funs_WKNSMSE.R` mse modules for running the category 1 SAM-based
    ICES MSY rule,
-   `MP_run.R` an R script for running MSE projections with the genetic
    algorithm, usually called from a job submission script (`*.pbs`),
-   `MP_analysis.R` script for analysing the results.

The following input files are provided:

-   `input/model_input_stk_d.rds` contains the XSA stock input data,
-   `input/model_input_idx.rds` contains the XSA survey index input
    data,
-   `input/stock.rds` contains the XSA model fit from WGCSE 2021 for the
    landings only assessment,
-   `input/stock_d.rds` contains the XSA model fit from WGCSE 2021 for
    the total catch assessment,
-   `input/fit.rds` SAM model fit,
-   `input/length/ALKs.csv` age-length keys, used for generating the
    length index,
-   `input/ple.27.7e/baseline/1000_20/input_rfb.rds` default mse input
    object for running the rfb rule.

## R, R packages and version info

The MSE simulations were run with R:

``` r
> sessionInfo()
R version 3.6.3 (2020-02-29)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19041)

Matrix products: default

locale:
[1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252   
[3] LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.1252    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggrepel_0.8.1          mse_2.0.3              FLBRP_2.5.4            data.table_1.12.8     
 [5] doParallel_1.0.15      tidyr_1.0.2            dplyr_0.8.4            foreach_1.4.8         
 [9] stockassessment_0.10.0 ggplotFL_2.6.8         FLfse_0.0.0.9007       FLash_2.5.11          
[13] FLXSA_2.6.4            FLAssess_2.6.3         FLCore_2.6.14.9004     iterators_1.0.13      
[17] lattice_0.20-40        ggplot2_3.3.3         

loaded via a namespace (and not attached):
 [1] xfun_0.12        tidyselect_1.0.0 TMB_1.7.19       purrr_0.3.3      splines_3.6.3    colorspace_2.0-0
 [7] vctrs_0.3.6      htmltools_0.4.0  stats4_3.6.3     yaml_2.2.1       mgcv_1.8-31      utf8_1.1.4      
[13] rlang_0.4.10     pillar_1.4.7     glue_1.3.1       withr_2.3.0      lifecycle_0.2.0  munsell_0.5.0   
[19] gtable_0.3.0     evaluate_0.14    codetools_0.2-16 knitr_1.28       labeling_0.4.2   fansi_0.4.1     
[25] Rcpp_1.0.5       scales_1.1.1     farver_2.0.3     gridExtra_2.3    ellipse_0.4.2    digest_0.6.27   
[31] stringi_1.4.6    grid_3.6.3       cli_2.2.0        tools_3.6.3      magrittr_2.0.1   tibble_2.1.3    
[37] crayon_1.3.4     pkgconfig_2.0.3  MASS_7.3-51.5    ellipsis_0.3.1   Matrix_1.2-18    rmarkdown_2.1   
[43] assertthat_0.2.1 rstudioapi_0.13  R6_2.5.0         nlme_3.1-144     compiler_3.6.3 
```

The framework is based on the Fisheries Library in R (FLR) framework.
The exact versions of the packages as used here can be installed with
`remotes`:

``` r
remotes::install_github(repo = "flr/FLCore", ref = "3d694903b9e6717b86c3e8486fc14ebf92908786")
remotes::install_github(repo = "shfischer/FLash", ref = "d1fb86fa081aaa5b6980d74b07d9adb44ad19a7f") # silenced version of FLash
remotes::install_github(repo = "flr/FLBRP", ref = "3a4d6390abc56870575fbaba3637091036468217")
remotes::install_github(repo = "flr/FLAssess", ref = "14c7752efcc5a2f463d0e9846be33b327c77b59b")
remotes::install_github(repo = "flr/FLXSA", ref = "40756e4d1eb60a0ab50d6e511477a8ee8dde2373")
remotes::install_github(repo = "flr/ggplotFL")
```

Furthermore, a data-limited fork of the `flr/mse` package is required:

``` r
remotes::install_github(repo = "shfischer/mse", ref = "mseDL2.0")
```

For the generation of the operating model, the SAM
[`stockassessment`](https://github.com/fishfollower/SAM/) R package from
GitHub is required:

-   `stockassessment` 0.10.0

Install this version with

``` r
install.packages("TMB") ### required to use stockassessment
devtools::install_github("fishfollower/SAM/stockassessment", ref = "6ddb31ffee0d867c569a86cd4dc7a9ba8243764c")
```

In order to use SAM within FLR, the functionality of the inofficial
package `FLfse` is used:

-   `FLfse` 1.0

``` r
devtools::install_github("shfischer/FLfse/FLfse", ref = "64cec8ee86c9309a5b4e84ff1790d0ed594bc7ee")
```

Furthermore, some more R packages available from CRAN are required:

``` r
install.packages(c("foreach", "DoParallel", "doRNG", "dplyr", "tidyr", "ggplot2", "Cairo", "scales")) 
```
