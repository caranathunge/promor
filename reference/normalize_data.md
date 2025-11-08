# Normalize intensity data

This function normalizes data using a user-specified normalization
method.

## Usage

``` r
normalize_data(df, method = "quantile")
```

## Arguments

- df:

  An `imp_df` object with missing values imputed using `impute_na` or a
  `raw_df` object containing missing values.

- method:

  Name of the normalization method to use. Choices are
  `"none", "scale", "quantile" or "cyclicloess."` Default is
  `"quantile."`

## Value

A `norm_df` object, which is a data frame of normalized protein
intensities.

## Details

- This function normalizes intensity values to achieve consistency among
  samples.

- It assumes that the intensities in the data frame have been
  log-transformed, therefore, it is important to make sure that
  `create_df` was run with `log_tr = TRUE`(default) when creating the
  `raw_df` object.

## See also

- `impute_na`

- See
  [`normalizeBetweenArrays`](https://rdrr.io/pkg/limma/man/normalizebetweenarrays.html)
  in the R package `limma` for more information on the different
  normalization methods available.

## Author

Chathurani Ranathunge

## Examples

``` r
## Generate a raw_df object with default settings. No technical replicates.
raw_df <- create_df(
  prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
  exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt"
)
#> 0 empty row(s) removed.
#> 0 empty column(s) removed.
#> 80 protein(s) (rows) only identified by site removed.
#> 65 reverse protein(s) (rows) removed.
#> 42 protein potential contaminant(s) (rows) removed.
#> 1923 protein(s) identified by 2 or fewer unique peptides removed.
#> Zeros have been replaced with NAs.
#> Data have been log-transformed.

## Impute missing values in the data frame using the default minProb
## method prioir to normalization.
imp_df <- impute_na(raw_df)

## Normalize the imp_df object using the default quantile method
norm_df1 <- normalize_data(imp_df)

## Use the cyclicloess method
norm_df2 <- normalize_data(imp_df, method = "cyclicloess")

## Normalize data in the raw_df object prior to imputation.
norm_df3 <- normalize_data(raw_df)
```
