
<!-- README.md is generated from README.Rmd. Please edit that file -->

# indoortransformer

<!-- badges: start -->
<!-- badges: end -->

The *indoortransformer* package is designed to extend the functionality
of existing commercial chemical databases by predicting likely indoor
transformation products. The current build is designed to be used only
for organophosphorus compounds (OPCs).

## Installation

The ChemmineR package is required to use indoortransformer. You can
install ChemmineR with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChemmineR")
```

The [latest version of
Rtools](https://cran.r-project.org/bin/windows/Rtools/) is also required
in order to install packages directly from Github.

You can install the current version of indoortransformer from GitHub
with:

``` r
install.packages("devtools")
devtools::install_github("toxicpeng/indoortransformer")
```

## Troubleshooting

If the package does not install correctly, first make sure that your
version of R is up to date (&gt;= 4.1.2).

If you encounter the following error while trying to install this
package, you may need to update/re-install Java on your computer.

``` r
# Error: .onLoad failed in loadNamespace() for 'rJava', details:
#   call: fun(libname, pkgname)
#   error: JAVA_HOME cannot be determined from the Registry
```

[Download the appropriate version of the Java SDK for your system
here](https://www.oracle.com/java/technologies/downloads/). Close and
reopen R. You will then need to set “JAVA\_HOME” to the directory on
your computer which contains the file “jvm.dll”.

``` r
# Example file path:
Sys.setenv(JAVA_HOME = 'C:\\Program Files\\Java\\jdk-17.0.1\\bin\\server')
```

If you are still unable to install indoortransformer, use the following
instead, then try to install indoortransformer again.

``` r
Sys.setenv(JAVA_HOME = '')
```

## Example

The *trans.Products()* function accepts a character vector of molecules
in SMILES (Simplified Molecular-Input Line-Entry System) format. For
example:

``` r
x <- c("CO[P](=S)(OC)Oc1ccc(cc1)[S](=O)(=O)N(C)C","CO[P](=O)(OC)OC=C(Cl)Cl")
y <- trans.Products(x)'
View(y)
```
