# Remove user-specified samples

This function removes user-specified samples from the data frame.

## Usage

``` r
rem_sample(raw_df, rem)
```

## Arguments

- raw_df:

  A `raw_df` object.

- rem:

  Name of the sample to remove.

## Value

A `raw_df` object.

## Details

- If all the technical replicates representing a sample needs to be
  removed, provide "Group_UniqueSampleID" as `rem`.

- If a specific technical replicate needs to be removed in case it shows
  weak correlation with other technical replicates for example, you can
  remove that particular technical replicate by providing
  "Group_UniqueSampleID_TechnicalReplicate" as `rem`.

## See also

[`corr_plot`](https://caranathunge.github.io/promor/reference/corr_plot.md),
[`create_df`](https://caranathunge.github.io/promor/reference/create_df.md)

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
# Check the first few rows of the raw_df object
head(raw_df)
#>                WT_4_1   WT_4_2   WT_4_3   WT_5_1   WT_5_2   WT_5_3   WT_6_1
#> gi|118496616 28.29364 28.23027 28.21688 28.06928 28.30029 28.16715 28.32976
#> gi|118496617 29.61863 29.73206 29.72970 29.43474 29.69021 29.79974 29.88464
#> gi|118496621       NA 23.96754       NA       NA       NA       NA       NA
#> gi|118496635 30.79452 30.80939 30.81917 30.59584 30.71367 30.74791 30.58924
#> gi|118496636 29.06789 29.18988 29.19237 28.99674 28.88547 28.99502 29.10345
#> gi|118496637 28.04062 27.74919 28.09763 27.91462 27.67364 27.69617 28.00544
#>                WT_6_2   WT_6_3   D8_1_1   D8_1_2   D8_1_3   D8_2_1   D8_2_2
#> gi|118496616 28.33360 28.30978 27.84458 28.18980 28.26794 28.22246 28.25953
#> gi|118496617 29.93382 29.96168 29.44843 29.89179 29.76792 29.85262 29.94182
#> gi|118496621       NA 24.31872       NA       NA       NA       NA       NA
#> gi|118496635 30.73036 30.64727 30.62366 30.64702 30.71857 30.53594 30.52739
#> gi|118496636 29.11018 28.95471 29.14163 28.98191 29.07324 28.94067 28.90216
#> gi|118496637 27.96560 27.96664 27.70352 27.72715 27.70227 27.80839 27.56005
#>                D8_2_3   D8_3_1   D8_3_2   D8_3_3
#> gi|118496616 28.34363 28.04004 28.30795 28.28639
#> gi|118496617 29.99010 30.20560 30.09723 30.21824
#> gi|118496621       NA       NA       NA       NA
#> gi|118496635 30.59290 30.66195 30.62793 30.60858
#> gi|118496636 28.85896 29.09419 28.89331 29.04336
#> gi|118496637 27.58992 27.63810 27.55866 27.46014

## Remove all technical replicates of "WT_4"
raw_df1 <- rem_sample(raw_df, "WT_4")

## Remove only technical replicate number 2 of "WT_4"
raw_df2 <- rem_sample(raw_df, "WT_4_2")
```
