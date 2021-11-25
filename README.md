
<!-- README.md is generated from README.Rmd. Please edit that file -->

# indoortransformer

<!-- badges: start -->
<!-- badges: end -->

The *indoortransformer* package is designed to extend the functionality
of existing commercial chemical databases by predicting likely indoor
transformation products. Build 0.9.0 is currently designed for use only
with organophosphorus compounds (OPCs).

## Installation

You can install the current version of indoortransformer from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("toxicpeng/indoortransformer")
```

## Example

The *trans.Products()* function accepts a character vector of molecules
in SMILES (Simplified Molecular-Input Line-Entry System) format. For
example:

``` r
# x <- c("CO[P](=S)(OC)Oc1ccc(cc1)[S](=O)(=O)N(C)C","CO[P](=O)(OC)OC=C(Cl)Cl")
# trans.Products(x)
```
