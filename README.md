
# immcp: Polypharmacology Toolkit for Traditional Chinese Medicine Research
[![](https://img.shields.io/cran/v/immcp?label=CRAN)](https://CRAN.R-project.org/package=immcp)
[![](https://img.shields.io/badge/devel%20version-1.0.4-blue)](https://github.com/YuanlongHu/immcp)
[![](https://img.shields.io/github/license/YuanlongHu/immcp)](https://github.com/YuanlongHu/immcp/blob/master/LICENSE.md)
[![](https://img.shields.io/github/workflow/status/YuanlongHu/immcp/R-CMD-check?label=macOS%2Fubuntu)](https://github.com/YuanlongHu/immcp/actions)
This R package was a toolkit for TCM polypharmacology research. Based on the biological descriptors and drug-disease interaction networks, it can analyze the potential polypharmacological mechanisms of TCM and be used for drug repositioning in TCM. 
+ Extract biological descriptors by defining the genesets, and calculate the similarity of the drug to the disease characterized by the biological descriptors
+ Create and analyze Drug-Disease Network


## Installation
You can install it from CRAN:

``` r
install.packages("immcp")
```

You can also install a development release from GitHub using the devtools package :

``` r
install.packages("devtools")
library(devtools)
install_github("YuanlongHu/immcp")
```

