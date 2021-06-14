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

## Quick-start

This is a basic using pre-processed data

``` r
library(coher)

data('coher_example')

result <- coher(coher_example$Y, coher_example$X)
```

## Example

This is a more thorough example includes the loading of variant and phenotypic data from a vcf and csv file.

Load libraries and data.

```r
library(coher)
library(vcfR)
library(ggplot2)

vcf.file.name <- system.file("extdata", "smalldat_coher.vcf.gz", package="coher")
vcf <- vcfR::read.vcfR(vcf.file.name)

photype.file.name <- system.file("extdata", "phenotypes.csv.gz", package="coher")
phenotypes <- read.csv(photype.file.name)
```

Convert data to format required by coher

```r
X <- vcfR::extract.gt(vcf)
#convert to binary
X <- t(1*(X=="1/1"))

Y <- lapply(2:ncol(phenotypes), function(j) setNames(phenotypes[, j], phenotypes$isolate))
```

Run coher

```r
result <- coher(Y, X)
result$coher.matrix
```

We can now plot the results using ggplot2

```r
plotdf <- data.frame(phenotype.A = rep(colnames(phenotypes)[2:ncol(phenotypes)], 3),
           phenotype.B = rep(colnames(phenotypes)[2:ncol(phenotypes)], each=3),
           coher = c(result$coher.matrix), stringsAsFactors = FALSE)
           
ggplot(plotdf, aes(x=phenotype.A, y=phenotype.B, fill=coher)) +
  geom_tile() +
  theme_bw(base_size = 14) +
  xlab('') + ylab('') +
  scale_fill_distiller(palette = 5, type='div')
```
