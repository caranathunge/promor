# Heatmap of differentially expressed proteins

This function generates a heatmap to visualize differentially expressed
proteins between groups

## Usage

``` r
heatmap_de(
  fit_df,
  df,
  adj_method = "BH",
  cutoff = 0.05,
  lfc = 1,
  sig = "adjP",
  n_top = 20,
  palette = "viridis",
  text_size = 10,
  save = FALSE,
  file_path = NULL,
  file_name = "HeatmapDE",
  file_type = "pdf",
  dpi = 80,
  plot_height = 7,
  plot_width = 7
)
```

## Arguments

- fit_df:

  A `fit_df` object from performing `find_dep`.

- df:

  The `norm_df` object or the `imp_df` object from which the `fit_df`
  object was obtained.

- adj_method:

  Method used for adjusting the p-values for multiple testing. Default
  is `"BH"`.

- cutoff:

  Cutoff value for p-values and adjusted p-values. Default is 0.05.

- lfc:

  Minimum absolute log2-fold change to use as threshold for differential
  expression. Default is 1.

- sig:

  Criteria to denote significance. Choices are `"adjP"` (default) for
  adjusted p-value or `"P"` for p-value.

- n_top:

  Number of top hits to include in the heat map.

- palette:

  Viridis color palette option for plots. Default is `"viridis"`. See
  [`viridis`](https://sjmgarnier.github.io/viridisLite/reference/viridis.html)
  for available options.

- text_size:

  Text size for axis text, labels etc.

- save:

  Logical. If `TRUE` saves a copy of the plot in the directory provided
  in `file_path`.

- file_path:

  A string containing the directory path to save the file.

- file_name:

  File name to save the plot. Default is "HeatmapDE."

- file_type:

  File type to save the plot. Default is `"pdf".`

- dpi:

  Plot resolution. Default is `80.`

- plot_height:

  Height of the plot. Default is 7.

- plot_width:

  Width of the plot. Default is 7.

## Value

A `ggplot2` plot object.

## Details

By default the tiles in the heatmap are reordered by intensity values
along both axes (x axis = samples, y axis = proteins).

## See also

- `find_dep`

- [`topTable`](https://rdrr.io/pkg/limma/man/toptable.html) and
  [`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html) functions from the
  [`limma`](https://rdrr.io/pkg/limma/man/01Introduction.html) package.

## Author

Chathurani Ranathunge

## Examples

``` r
## Build a heatmap of differentially expressed proteins using the provided
## example fit_df and norm_df data objects
heatmap_de(covid_fit_df, covid_norm_df)
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the promor package.
#>   Please report the issue at <https://github.com/caranathunge/promor/issues>.


## Create a heatmap with P-value of 0.05 and log fold change of 1 as
## significance criteria.
heatmap_de(covid_fit_df, covid_norm_df, cutoff = 0.05, sig = "P")


## Visualize the top 30 differentially expressed proteins in the heatmap and
## change the color palette
heatmap_de(covid_fit_df, covid_norm_df,
  cutoff = 0.05, sig = "P", n_top = 30,
  palette = "magma"
)

```
