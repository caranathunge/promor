# Identify differentially expressed proteins between groups

This function performs differential expression analysis on protein
intensity data with limma.

## Usage

``` r
find_dep(
  df,
  save_output = FALSE,
  save_tophits = FALSE,
  file_path = NULL,
  adj_method = "BH",
  cutoff = 0.05,
  lfc = 1,
  n_top = 20
)
```

## Arguments

- df:

  A `norm_df` object or an `imp_df` object.

- save_output:

  Logical. If `TRUE` saves results from the differential expression
  analysis in a text file labeled "limma_output.txt" in the directory
  specified by `file_path`.

- save_tophits:

  Logical. If `TRUE` saves `n_top` number of top hits from the
  differential expression analysis in a text file labeled "TopHits.txt"
  in the directory specified by `file_path`.

- file_path:

  A string containing the directory path to save the file.

- adj_method:

  Method used for adjusting the p-values for multiple testing. Default
  is `"BH"` for "Benjamini-Hochberg" method.

- cutoff:

  Cutoff value for p-values and adjusted p-values. Default is 0.05.

- lfc:

  Minimum absolute log2-fold change to use as threshold for differential
  expression.

- n_top:

  The number of top differentially expressed proteins to save in the
  "TopHits.txt" file. Default is `20`.

## Value

A `fit_df` object, which is similar to a `limma` `fit` object.

## Details

- `save_output` saves the complete results table from the differential
  expression analysis.

- `save_tophits` first subsets the results to those with absolute log
  fold change of more than 1, performs multiple correction with the
  method specified in `adj_method` and outputs the top `n_top` results
  based on lowest p-value and adjusted p-value.

- If the number of hits with absolute log fold change of more than 1 is
  less than `n_top`, `find_dep` prints only those with log-fold change
  \> 1 to "TopHits.txt".

- If the `file_path` is not specified, text files will be saved in a
  temporary directory.

## References

Ritchie, Matthew E., et al. "limma powers differential expression
analyses for RNA-sequencing and microarray studies." Nucleic acids
research 43.7 (2015): e47-e47.

## See also

- [`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html),
  [`eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html),
  [`topTable`](https://rdrr.io/pkg/limma/man/toptable.html), and
  [`write.fit`](https://rdrr.io/pkg/limma/man/writefit.html) functions
  from the [`limma`](https://rdrr.io/pkg/limma/man/01Introduction.html)
  package.

## Author

Chathurani Ranathunge

## Examples

``` r
## Perform differential expression analysis using default settings
fit_df1 <- find_dep(ecoli_norm_df)
#> 1186 siginificantly differentially expressed proteins found.

## Change p-value and adjusted p-value cutoff
fit_df2 <- find_dep(ecoli_norm_df, cutoff = 0.1)
#> 1227 siginificantly differentially expressed proteins found.
```
