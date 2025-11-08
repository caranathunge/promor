# Pre-process protein intensity data for modeling

This function pre-processes protein intensity data from the top
differentially expressed proteins identified with `find_dep` for
modeling.

## Usage

``` r
pre_process(
  fit_df,
  norm_df,
  sig = "adjP",
  sig_cutoff = 0.05,
  fc = 1,
  n_top = 20,
  find_highcorr = TRUE,
  corr_cutoff = 0.9,
  save_corrmatrix = FALSE,
  file_path = NULL,
  rem_highcorr = TRUE
)
```

## Arguments

- fit_df:

  A `fit_df` object from performing `find_dep`.

- norm_df:

  The `norm_df` object from which the `fit_df` object was obtained.

- sig:

  Criteria to denote significance in differential expression. Choices
  are `"adjP"` (default) for adjusted p-value or `"P"` for p-value.

- sig_cutoff:

  Cutoff value for p-values and adjusted p-values in differential
  expression. Default is `0.05`.

- fc:

  Minimum absolute log-fold change to use as threshold for differential
  expression. Default is `1`.

- n_top:

  The number of top hits from `find_dep` to be used in modeling. Default
  is `20`.

- find_highcorr:

  Logical. If `TRUE` (default), finds highly correlated proteins.

- corr_cutoff:

  A numeric value specifying the correlation cutoff. Default is `0.90`.

- save_corrmatrix:

  Logical. If `TRUE`, saves a copy of the protein correlation matrix in
  a tab-delimited text file labeled "Protein_correlation.txt" in the
  directory specified by `file_path`.

- file_path:

  A string containing the directory path to save the file.

- rem_highcorr:

  Logical. If `TRUE` (default), removes highly correlated proteins
  (predictors or features).

## Value

A `model_df` object, which is a data frame of protein intensities with
proteins indicated by columns.

## Details

This function creates a data frame that contains protein intensities for
a user-specified number of top differentially expressed proteins.

- Note: Most models will benefit from reducing correlation between
  proteins (predictors or features), therefore we recommend removing
  those proteins at this stage to reduce pairwise-correlation.

- If no or few proteins meet the significance threshold for differential
  expression, you may adjust `sig`, `fc`, and/or `sig_cutoff`
  accordingly to obtain more proteins for modeling.

## See also

- `find_dep`, `normalize_data`

- [`caret: findCorrelation`](https://rdrr.io/pkg/caret/man/findCorrelation.html)

## Author

Chathurani Ranathunge

## Examples

``` r
## Create a model_df object with default settings.
covid_model_df1 <- pre_process(fit_df = covid_fit_df, norm_df = covid_norm_df)
#> Total number of differentially expressed proteins (8) is less than n_top.
#> None of the proteins show high pair-wise correlation.
#> 
#> No highly correlated proteins to be removed.

## Change the correlation cutoff.
covid_model_df2 <- pre_process(covid_fit_df, covid_norm_df, corr_cutoff = 0.95)
#> Total number of differentially expressed proteins (8) is less than n_top.
#> None of the proteins show high pair-wise correlation.
#> 
#> No highly correlated proteins to be removed.

## Change the significance criteria to include more proteins
covid_model_df3 <- pre_process(covid_fit_df, covid_norm_df, sig = "P")
#> None of the proteins show high pair-wise correlation.
#> 
#> No highly correlated proteins to be removed.

## Change the number of top differentially expressed proteins to include
covid_model_df4 <- pre_process(covid_fit_df, covid_norm_df, sig = "P", n_top = 24)
#> None of the proteins show high pair-wise correlation.
#> 
#> No highly correlated proteins to be removed.
```
