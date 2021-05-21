# coher

<!-- badges: start -->
<!-- badges: end -->

Coher estimates the genetic correlations between different phenotypes in bacterial genome association studies.

## Installation

You can install the development version of coher by running

``` r
install.packages('remotes')
remotes::install_github("tienmt/coher")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
data('coher_example')
result <- coher(coher_example$Y, coher_example$X)
```
