
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SimTOST <img src="man/figures/logo.png" align="right" width="150" alt="" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/smartdata-analysis-and-statistics/SimTOST/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/smartdata-analysis-and-statistics/SimTOST/actions/workflows/R-CMD-check.yml)
[![CRAN
status](https://www.r-pkg.org/badges/version/SimTOST)](https://CRAN.R-project.org/package=SimTOST)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/last-month/SimTOST)](https://cran.r-project.org/package=SimTOST)
<!-- badges: end -->

`SimTOST` is an R package specifically designed for bioequivalence
studies, providing simulation-based sample size estimation for the Two
One-Sided Tests (TOST) procedure. It offers flexible options to handle
complex study designs, including trials with multiple correlated primary
endpoints, multiple hypotheses, and treatment arms. By incorporating
correlations between endpoints, `SimTOST` ensures accurate and robust
planning of bioequivalence trials, making it a powerful tool for studies
with intricate requirements.

## Installation

`SimTOST` is available on CRAN and can be installed by running the
following code.

``` r
install.packages("SimTOST")
```

You can also install the development version of SimTOST from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("smartdata-analysis-and-statistics/SimTOST")
```

## Vignettes

The main features of this package is `sampleSize` function which can be
used to calculate sample size for individual and multiple endpoints.
Various worked examples are available as
[vignettes](https://smartdata-analysis-and-statistics.github.io/SimTOST/),
with an introduction provided at
[Introduction](https://smartdata-analysis-and-statistics.github.io/SimTOST/articles/intropkg.html).

| **Vignette** | **Design** | **Number of Arms** | **Number of Endpoints** |
|----|----|----|----|
| [Bioequivalence Tests for Parallel Trial Designs with Log-Normal Data](https://smartdata-analysis-and-statistics.github.io/SimTOST/articles/sampleSize_parallel.html) | Parallel | 2 | Multiple (e.g., 2 or 3) |
| [Bioequivalence Tests for 2x2 Cross-Over Trial Designs with Log-Normal Data](https://smartdata-analysis-and-statistics.github.io/SimTOST/articles/sampleSize_crossover.html) | Cross-over | 2 | 2 |
| [Bioequivalence Tests for Parallel Trial Designs: 2 Arms, 1 Endpoint](https://smartdata-analysis-and-statistics.github.io/SimTOST/articles/sampleSize_parallel_2A1E.html) | Parallel | 2 | 1 |
| [Bioequivalence Tests for Parallel Trial Designs: 3 Arms, 1 Endpoint](https://smartdata-analysis-and-statistics.github.io/SimTOST/articles/sampleSize_parallel_3A1E.html) | Parallel | 3 | 1 |
| [Bioequivalence Tests for Parallel Trial Designs: 2 Arms, 3 Endpoints](https://smartdata-analysis-and-statistics.github.io/SimTOST/articles/sampleSize_parallel_2A3E.html) | Parallel | 2 | 3 |
| [Bioequivalence Tests for Parallel Trial Designs: 3 Arms, 3 Endpoints](https://smartdata-analysis-and-statistics.github.io/SimTOST/articles/sampleSize_parallel_3A3E.html) | Parallel | 3 | 3 |
