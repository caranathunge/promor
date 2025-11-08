# Remove user-specified proteins (features) from a data frame

This function removes user-specified proteins from a `model_df` object

## Usage

``` r
rem_feature(model_df, rem_protein)
```

## Arguments

- model_df:

  A `model_df` object.

- rem_protein:

  Name of the protein to remove.

## Value

A `model_df` object.

## Details

- For example, you can choose to remove a protein from the `model_df`
  object if the protein does not show distinct patterns of variation
  among conditions. This protein may show mostly overlapping
  distributions in the feature plots.

- Another incidence would be removing a protein that is very low in
  variable importance in the models built using `train_models`. You can
  visualize variable importance using `varimp_plot`.

## See also

[`feature_plot`](https://caranathunge.github.io/promor/reference/feature_plot.md),
[`pre_process`](https://caranathunge.github.io/promor/reference/pre_process.md)

## Author

Chathurani Ranathunge

## Examples

``` r
covid_model_df <- pre_process(fit_df = covid_fit_df, norm_df = covid_norm_df)
#> Total number of differentially expressed proteins (8) is less than n_top.
#> None of the proteins show high pair-wise correlation.
#> 
#> No highly correlated proteins to be removed.

## Remove sp|P22352|GPX3_HUMAN protein from the model_df object
covid_model_df1 <- rem_feature(covid_model_df, rem_protein = "sp|P22352|GPX3_HUMAN")
#> Protein sp|P22352|GPX3_HUMAN has been removed.
```
