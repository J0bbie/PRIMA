# PRIMA

![license](https://img.shields.io/badge/license-GPL--3-blue.svg) [![GitHub issues](https://img.shields.io/github/issues/J0bbie/PRIMA.svg)]() ![rversion](https://img.shields.io/badge/R%20version-%3E4.1.2-lightgrey.svg)

# Introduction

`PRIMA` is a wrapper-package around the use of [primer3](https://github.com/primer3-org/primer3) for designing PCR primers.
It contains several functions to output only those primer3-designed primers which include the given target DNA sequence.

This is useful to ensure that the resulting PCR-product contains for instance the somatic breakpoint of structural variants or other necessary genomic variants-of-interest.


# Installation

The latest development version can be installed directly from GitHub:

```R
# Require/install devtools package if not already installed.
if (!require("devtools")) install.packages("devtools", repos = "http://cran.r-project.org")

# Install ProteoDisco from GitHub.
devtools::install_github(repo = "J0bbie/PRIMA")

# Download all required packages
library(PRIMA)
```

# Usage

Please view the vignettes for instructions and tutorials on how to use this package.

```R
browseVignettes("PRIMA")
```
