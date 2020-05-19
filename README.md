
# immcp

<!-- badges: start -->
<!-- badges: end -->



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
data("drugSample")
FP <- extrFP(disease_biomarker = drugSample$disease_biomarker,
               drug_target = drugSample$herb_target,
               geneset = "ImmGenTop150")
res <- score_fp(FP, n=100)
res <- as.data.frame(res)
```

