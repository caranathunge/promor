# promor
### A comprehensive R package for proteomics data analysis and modeling

*promor* is a comprehensive R package that combines proteomics data analysis 
with machine learning approaches to assess the predictive performance of
specific proteins. The input for *promor* is the proteinGroups.txt file 
produced by *MaxQuant*. *promor* provides a wealth of data analysis and 
visualization tools at the protein level to analyze label-free proteomics data.
Further, it provides the user with
the option to test the predictive performance of proteins identified by
differential expression analysis using machine learning approaches.


## Installation

*promor*  can be installed directly from GitHub.

```r
install.packages("promor", dependencies = TRUE)

```

## Tutorials

We have provided several tutorial vignettes to help get you started with 
*promor*. 
We suggest starting with our 
to get a basic understanding of *promor* and its functions.
```r
vignette("promor", package="promor")
```
If your data contains technical replicates, you can refer to 
for specific *promor* functions that you can use in your analysis.
```r
vignette("promor", package="promor")
```

If your data does NOT contain technical replicates, you may follow the 
workflow provided in 

```r
vignette("promor", package="promor")
```


copies of the vignettes can be found here:

- [promor](http://promor.html)

## Citation
[Journal](https://www.blabla)
