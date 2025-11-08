# Compute average intensity

This function computes average intensities across technical replicates
for each sample.

## Usage

``` r
aver_techreps(raw_df)
```

## Arguments

- raw_df:

  A `raw_df` object containing technical replicates.

## Value

A `raw_df` object of averaged intensities.

## Details

`aver_techreps` assumes that column names in the data frame follow the
"Group_UniqueSampleID_TechnicalReplicate" notation. (Use `head(raw_df)`
to see the structure of the `raw_df` object.)

## See also

- [`create_df`](https://caranathunge.github.io/promor/reference/create_df.md)

## Author

Chathurani Ranathunge

## Examples

``` r
## Use a data set containing technical replicates to create a raw_df object
raw_df <- create_df(
prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg2.txt",
exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed2.txt",
tech_reps = TRUE
)
#> 0 empty row(s) removed.
#> 1 empty column(s) removed.
#> 12 reverse protein(s) (rows) removed.
#> 29 protein contaminant(s) (rows) removed.
#> 188 protein(s) identified by 2 or fewer unique peptides removed.
#> Zeros have been replaced with NAs.
#> Data have been log-transformed.

# Compute average intensities across technical replicates.
rawdf_ave <- aver_techreps(raw_df)
```
