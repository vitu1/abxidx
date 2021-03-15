
<!-- badges: start -->

[![R-CMD-check](https://github.com/vitu1/abxidx/workflows/R-CMD-check/badge.svg)](https://github.com/vitu1/abxidx/actions)
[![Travis build
status](https://travis-ci.com/vitu1/abxidx.svg?branch=main)](https://travis-ci.com/vitu1/abxidx)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/vitu1/abxidx?branch=master&svg=true)](https://ci.appveyor.com/project/VincentTu/abxidx)

<!-- badges: end -->

# abxidx

The goal of abxidx is to calculate an index for a given bacterial
community’s susceptibility to the specified
antibiotics

<!-- ## Installation -->

<!-- You can install the released version of abxidx from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("abxidx") -->

<!-- ``` -->

## Example

``` r
library(abxidx)

vancomycin_index(abx_test_df)
#>         a         b         c         d 
#> 0.0000000 0.4771213      -Inf       Inf
```
