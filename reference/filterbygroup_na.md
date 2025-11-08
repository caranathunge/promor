# Filter proteins by group level missing data

This function filters out proteins based on missing data at the group
level.

## Usage

``` r
filterbygroup_na(raw_df, set_na = 0.34, filter_condition = "either")
```

## Arguments

- raw_df:

  A `raw_df` object (output of
  [`create_df`](https://caranathunge.github.io/promor/reference/create_df.md))

- set_na:

  The proportion of missing data allowed. Default is 0.34 (one third of
  the samples in the group).

- filter_condition:

  If set to `"each"`, proteins that exceed the missing value proportion
  threshold set by `set_na` in each group will be removed (lenient). If
  set to `"either"`(default), proteins that exceed the missing value
  proportion threshold set by `set_na` in at least one group will be
  removed (stringent).

## Value

A `raw_df` object.

## Details

- If `filter_condition = "each"`, it then removes proteins (rows) from
  the data frame if the proportion of NAs in **each** group exceeds the
  threshold indicated by `set_na` (default is 0.34). This option is more
  lenient in comparison to `filter_condition = "either"`, where proteins
  that exceeds the missing data threshold in **either** group gets
  removed from the data frame.

## See also

[`create_df`](https://caranathunge.github.io/promor/reference/create_df.md)

## Author

Chathurani Ranathunge

## Examples

``` r
# \donttest{
# Generate a raw_df object with default settings. No technical replicates.
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

## Remove proteins that exceed 34% NAs in either group (default)
rawdf_filt1 <- filterbygroup_na(raw_df)
#> 224 proteins with higher than 34% NAs in at least one group removed.

## Remove proteins that exceed 34% NAs in each group
rawdf_filt2 <- filterbygroup_na(raw_df, filter_condition = "each")
#> 65 proteins with higher than 34% NAs in each group removed.

## Proportion of samples with NAs allowed in each group = 0.5
rawdf_filt3 <- filterbygroup_na(raw_df, set_na = 0.5, filter_condition = "each")
#> 65 proteins with higher than 50% NAs in each group removed.

# }
```
