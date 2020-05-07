
# immcp

<!-- badges: start -->
<!-- badges: end -->

The goal of immcp is to discovery candidate prescriptions.

## Installation

You can install it from GitHub using the devtools package :

``` r
install.packages("devtools")
library(devtools)
install_github("YuanlongHu/immcp")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(immcp)
# demo data
data("COVID19")
# basic example code
res <- (disease_biomarker = COVID19$disease_biomarker,
        disease_network = COVID19$disease_network,
        target = COVID19$herb_target,
        geneset = NULL)
# adjust score
res_adj <- score_adjust(result = res, n = 100)
```

