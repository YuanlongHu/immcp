
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
res <- score_immcp()
# adjust score
res_adj <- score_adjust()
```

