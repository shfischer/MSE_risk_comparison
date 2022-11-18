Risk equivalence in data-limited and data-rich fisheries management: an
example based on the ICES advice framework
================

## Introduction

This repository contains the data and code for a comparison of the
data-rich ICES MSY rule (category 1) to the ICES category 3 data-limited
empirical methods (rfb rule, hr rule,
<https://doi.org/10.17895/ices.advice.19801564>) and is the basis for
the publication:

> Fischer, S. H., De Oliveira, J. A. A., Mumford, J. D., and Kell, L. T.
> (accepted). Risk equivalence in data-limited and data-rich fisheries
> management: an example based on the ICES advice framework. Fish and
> Fisheries.

Three case study stocks are included:

1.  Plaice (*Pleuronectes platessa*) in Division 7.e (western English
    Channel) (<https://doi.org/10.17895/ices.advice.7822>)

2.  Cod (*Gadus morhua*) in Subarea 4, Division 7.d, and Subdivision 20
    (North Sea, eastern English Channel, Skagerrak)
    (<https://doi.org/10.17895/ices.advice.7746>)

3.  Herring (*Clupea harengus*) in Subarea 4 and divisions 3.a and 7.d,
    autumn spawners (North Sea, Skagerrak and Kattegat, eastern English
    Channel) (<https://doi.org/10.17895/ices.advice.7770>)

The operating models (OMs) are created using the SAM
[`stockassessment`](https://github.com/fishfollower/SAM/) R package and
follows the principles developed during the ICES Workshop on North Sea
stocks management strategy evaluation
([WKNSMSE](https://doi.org/10.17895/ices.pub.5090)).

The simulation is based on the Fisheries Library in R
([FLR](http://www.flr-project.org/)) and its
[`mse`](https://github.com/flr/mse) package.

The data-limited MSE framework is an adaptation of Fischer et
al. ([2021](https://dx.doi.org/10.1093/icesjms/fsab169)), see the
[`shfischer/GA_MSE_PA`](https://github.com/shfischer/GA_MSE_PA) GitHub
repository for the original source code and documentation.

For exact reproducibility, R packages versions are recorded with
[renv](https://rstudio.github.io/renv/) in a
[renv.lock](https://github.com/shfischer/MSE_risk_comparison/blob/master/renv.lock)
file.

## Repository structure

-   `funs_*`: Function libraries, defining the functions used in the
    other scripts

    -   [`funs.R`](https://github.com/shfischer/MSE_risk_comparison/blob/master/funs.R):
        generic function library, including definition of data-limited
        management procedures (MPs)

    -   [`funs_analysis.R`](https://github.com/shfischer/MSE_risk_comparison/blob/master/funs_analysis.R):
        for analysis of results

    -   [`funs_GA.R`](https://github.com/shfischer/MSE_risk_comparison/blob/master/funs_GA.R):
        functions used in the optimisation with the genetic algorithm
        (GA)

    -   [`funs_OM.R`](https://github.com/shfischer/MSE_risk_comparison/blob/master/funs_OM.R):
        functions for generating the operating models

    -   [`funs_WKNSMSE.R`](https://github.com/shfischer/MSE_risk_comparison/blob/master/funs_WKNSMSE.R):
        functions required for the ICES MSY rule

-   `OM_*`: Scripts for operating models (OMs, including alternative
    OMs)

    -   [`OM_ple.27.7e.R`](https://github.com/shfischer/MSE_risk_comparison/blob/master/OM_ple.27.7e.R)
        for plaice

    -   [`OM_cod.27.47d20.R`](https://github.com/shfischer/MSE_risk_comparison/blob/master/OM_cod.27.47d20.R)
        for cod

    -   [`OM_her.27.3a47d.R`](https://github.com/shfischer/MSE_risk_comparison/blob/master/OM_her.27.3a47d.R)
        for herring

    -   [`OM_MSY.R`](https://github.com/shfischer/MSE_risk_comparison/blob/master/OM_MSY.R):
        script for estimating MSY

    -   [`OM_MSY.pbs`](https://github.com/shfischer/MSE_risk_comparison/blob/master/OM_MSY.pbs):
        job script, calling `OM_MSY.R`

-   `MP_*`: Script for running and analysing the MSE

    -   [`MP_analysis.R`](https://github.com/shfischer/MSE_risk_comparison/blob/master/MP_analysis.R):
        script for analysing MSE results (summarising, exporting,
        visualisation, …)

    -   [`MP_run.R`](https://github.com/shfischer/MSE_risk_comparison/blob/master/MP_run.R):
        script for running any MP in the MSE and optimising MPs

    -   `MP_*.pbs`: job submission scripts, used for running MP_run.R on
        a high-performance computing cluster,
        e.g. [`MP_run_rfb_mult.pbs`](https://github.com/shfischer/MSE_risk_comparison/blob/master/MP_run_rfb_mult.pbs)
        for optimising the rfb rule with a multiplier

-   [`MP_check_SPiCT.R`](https://github.com/shfischer/MSE_risk_comparison/blob/master/MP_check_SPiCT.R):
    script for checking the [SPiCT](https://github.com/DTUAqua/spict)
    model

[`input/`](https://github.com/shfischer/MSE_risk_comparison/tree/master/input):
This directory contains all files required for generating the OMs for
the three stocks (`OM_*.R`)

-   [`input/ple.27.7e/preparation/`](https://github.com/shfischer/MSE_risk_comparison/tree/master/input/ple.27.7e/preparation):
    for plaice

-   [`input/cod.27.47d20/preparation/`](https://github.com/shfischer/MSE_risk_comparison/tree/master/input/cod.27.47d20/preparation):
    for cod

-   [`input/her.27.3a47d/preparation/`](https://github.com/shfischer/MSE_risk_comparison/tree/master/input/her.27.3a47d/preparation):
    for herring

-   [`input/OM_refpts.csv`](https://github.com/shfischer/MSE_risk_comparison/blob/master/input/OM_refpts.csv):
    summarised OM reference points

[`output/`](https://github.com/shfischer/MSE_risk_comparison/tree/master/output):
This directory contains some summarised results

-   [`output/MPs_baseline.csv`](https://github.com/shfischer/MSE_risk_comparison/blob/master/output/MPs_baseline.csv):
    summary statistics and parameterisations of default and optimised
    MPs for all stocks

-   [`output/plots/wormplots/`](https://github.com/shfischer/MSE_risk_comparison/tree/master/output/plots/wormplots):
    Projections (wormplots) for all default and optimised MPs, for all
    stocks and OMs

## R, R packages and version info

The MSE simulations were run with R:

``` r
> sessionInfo()
R version 4.1.0 (2021-05-18)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19041)
```

The package versions and their dependencies are recorded with the R
package [renv](https://rstudio.github.io/renv/) and stored in the file
[renv.lock](https://github.com/shfischer/MSE_risk_comparison/blob/master/renv.lock).
The exact package version can be restored by cloning this repository,
navigating into this folder in R (or setting up a project), installing
the renv package

``` r
install.packages("renv")
```

and calling

``` r
renv::restore()
```

See [renv](https://rstudio.github.io/renv/) and the package
documentation for details.

The framework is based on the Fisheries Library in R (FLR) framework and
uses the [FLR packages](https://flr-project.org/)
[`FLCore`,](https://github.com/flr/FLCore)
[`FLasher`](https://github.com/flr/FLasher),
[`FLBRP`](https://github.com/flr/FLBRP),
[`FLAssess`](https://github.com/flr/FLAssess),
[`FLXSA`](https://github.com/flr/FLXSA),
[`ggplotFL`](https://github.com/flr/ggplotFL),
[`mse`](https://github.com/flr/mse), and
[`FLfse`](https://github.com/shfischer/FLfse). See
[renv.lock](https://github.com/shfischer/MSE_risk_comparison/blob/master/renv.lock)
for version details and sources.

Also, the R package
[`stockassessment`](https://github.com/fishfollower/SAM)is used.

For running the optimisations on a high-performance computing cluster, a
suitable MPI back-end and the R package
[`Rmpi`](https://cran.r-project.org/web/packages/Rmpi/index.html) are
required.
