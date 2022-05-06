
-   [promor](#promor)
    -   [Installation](#installation)
    -   [Proteomics data analysis with
        promor](#proteomics-data-analysis-with-promor)
    -   [Tutorials](#tutorials)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# promor

### A comprehensive R package for proteomics data analysis and modeling

<!-- badges: start -->
<!-- badges: end -->

*promor* is an R package that combines proteomics data analysis with
machine learning-based modeling. Input files for *promor* are the
*proteinGroups.txt* file produced by *MaxQuant* and an *expDesign.txt*
file, which contains the experimental design of your proteomics data.  
*promor* provides a wealth of data analysis and visualization tools at
the protein level to analyze label-free proteomics data.

## Installation

You can install the development version of promor from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caranathunge/promor")
```

Alternatively, *promor* can be installed directly from a source. First,
you need to download and save the package files in a local directory.
Then you may install *promor* as follows:

``` r
install.packages("path/to/promor", repos = NULL, type = "source", dependencies = TRUE)
```

Note: Full path and file name should be provided as “path/to/promor.”
For example, if promor is saved on your C drive, it could be :
`"C://promor_0.1.0.tar.gz"`

## Proteomics data analysis with promor

![promor flowchart by caranathunge](./promor_ProtAnalysis_BuYe.png)
*Figure 1. A schematic diagram of suggested workflows for proteomics
data analysis with promor.*

### Example

Here is a minimal working example that shows how to identify
differentially expressed proteins between two conditions in five simple
steps using *promor*. We use a previously published data set from [Cox
et al. (2014)](https://europepmc.org/article/MED/24942700#id609082).

``` r
#Load promor
library(promor)
#Create a raw.df object with the files provided in extdata folder.
raw <- create_df(prot.groups = system.file("extdata", "ecoli_proteinGroups.txt", 
                                           package = "promor"), 
                 exp.design = system.file("extdata", "expDesign.txt", 
                                          package = "promor"))
#Filter out proteins with high levels of missing data in each condition
raw_filtered <- filterbygroup_na(raw)
#Impute missing data and create an imp.df object.
imp_df <- impute_na(raw_filtered)
#Normalize data and create a norm.df object
norm_df <- normalize_data(imp_df)
#Perform differential expression analysis and create a fit.df object
fit_df <- find_dep(norm_df)
```

Lets take a look at the results using a volcano plot.

``` r
 volcano_plot(fit_df)
```

<img src="man/figures/README-volcanoplot-1.png" width="90%" style="display: block; margin: auto;" />

## Tutorials

You can choose a tutorial from the list below that best fits your
experiment and the structure of your proteomics data.

1.  This README file can be accessed from inside *promor* as follows,

``` r
vignette("intro_to_promor", package = "promor")
```

2.  If your data does NOT contain technical replicates, you can refer to
    the following tutorial.

``` r
vignette("promor_with_notechreps", package = "promor")
```

3.  If your data contains technical replicates, you can refer to the
    following tutorial for an illustrative example.

``` r
vignette("promor_with_techreps", package = "promor")
```
