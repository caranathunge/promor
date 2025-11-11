# Impute missing values

This function imputes missing values using a user-specified imputation
method.

## Usage

``` r
impute_na(
  df,
  method = "minProb",
  tune_sigma = 1,
  q = 0.01,
  maxiter = 10,
  ntree = 20,
  n_pcs = 2,
  seed = NULL
)
```

## Arguments

- df:

  A `raw_df` object (output of
  [`create_df`](https://caranathunge.github.io/promor/reference/create_df.md))
  containing missing values or a `norm_df` object after performing
  normalization.

- method:

  Imputation method to use. Default is `"minProb"`. Available methods:
  `"minDet", "RF", "kNN", and "SVD"`.

- tune_sigma:

  A scalar used in the `"minProb"` method for controlling the standard
  deviation of the Gaussian distribution from which random values are
  drawn for imputation.  
  Default is 1.

- q:

  A scalar used in `"minProb"` and `"minDet"` methods to obtain a low
  intensity value for imputation. `q` should be set to a very low value.
  Default is 0.01.

- maxiter:

  Maximum number of iterations to be performed when using the `"RF"`
  method. Default is `10`.

- ntree:

  Number of trees to grow in each forest when using the `"RF"` method.
  Default is `20`.

- n_pcs:

  Number of principal components to calculate when using the `"SVD"`
  method. Default is 2.

- seed:

  Numerical. Random number seed. Default is `NULL`

## Value

An `imp_df` object, which is a data frame of protein intensities with no
missing values.

## Details

- `impute_na` function imputes missing values using a user-specified
  imputation method from the available options, `minProb`, `minDet`,
  `kNN`, `RF`, and `SVD`.

- **Note: Some imputation methods may require that the data be
  normalized prior to imputation.**

- Make sure to fix the random number seed with `seed` for
  reproducibility

.

## References

Lazar, Cosmin, et al. "Accounting for the multiple natures of missing
values in label-free quantitative proteomics data sets to compare
imputation strategies." Journal of proteome research 15.4 (2016):
1116-1125.

## See also

More information on the available imputation methods can be found in
their respective packages.

- For `minProb` and `minDet` methods, see `imputeLCMD` package.

- For Random Forest (`RF`) method, see
  [`missForest`](https://rdrr.io/pkg/missForest/man/missForest.html).

- For `kNN` method, see [`kNN`](https://rdrr.io/pkg/VIM/man/kNN.html)
  from the [`VIM`](https://rdrr.io/pkg/VIM/man/VIM-package.html)
  package.

- For `SVD` method, see
  [`pca`](https://rdrr.io/pkg/pcaMethods/man/pca.html) from the
  [`pcaMethods`](https://rdrr.io/pkg/pcaMethods/man/pcaMethods.html)
  package.

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
## method.
imp_df1 <- impute_na(raw_df, seed = 3312)

# \donttest{
## Impute using the RF method with the number of iterations set at 5
## and number of trees set at 100.
imp_df2 <- impute_na(raw_df,
  method = "RF",
  maxiter = 5, ntree = 100,
  seed = 3312
)
#>   missForest iteration 1 in progress...done!
#>     estimated error(s): 0.1361906 
#>     difference(s): 0.0009932238 
#>     time: 10.544 seconds
#> 
#>   missForest iteration 2 in progress...done!
#>     estimated error(s): 0.1318301 
#>     difference(s): 2.294329e-06 
#>     time: 10.198 seconds
#> 
#>   missForest iteration 3 in progress...done!
#>     estimated error(s): 0.1318273 
#>     difference(s): 8.586608e-07 
#>     time: 9.886 seconds
#> 
#>   missForest iteration 4 in progress...done!
#>     estimated error(s): 0.1316124 
#>     difference(s): 7.427129e-07 
#>     time: 10.325 seconds
#> 
#>   missForest iteration 5 in progress...done!
#>     estimated error(s): 0.1314333 
#>     difference(s): 6.803003e-07 
#>     time: 10.111 seconds
#> 


## Using the kNN method.
imp_df3 <- impute_na(raw_df, method = "kNN", seed = 3312)
#>      H_2      H_3      L_1      L_2      L_3      H_2      H_3      L_1 
#> 18.85805 19.58752 19.10595 19.55137 18.57692 36.98375 37.25068 36.49326 
#>      L_2      L_3 
#> 36.27569 36.35529 
#>      H_1      H_3      L_1      L_2      L_3      H_1      H_3      L_1 
#> 18.48731 19.58752 19.10595 19.55137 18.57692 37.19013 37.25068 36.49326 
#>      L_2      L_3 
#> 36.27569 36.35529 
#>      H_1      H_2      L_1      L_2      L_3      H_1      H_2      L_1 
#> 18.48731 18.85805 19.10595 19.55137 18.57692 37.19013 36.98375 36.49326 
#>      L_2      L_3 
#> 36.27569 36.35529 
#>      H_1      H_2      H_3      L_2      L_3      H_1      H_2      H_3 
#> 18.48731 18.85805 19.58752 19.55137 18.57692 37.19013 36.98375 37.25068 
#>      L_2      L_3 
#> 36.27569 36.35529 
#>      H_1      H_2      H_3      L_1      L_3      H_1      H_2      H_3 
#> 18.48731 18.85805 19.58752 19.10595 18.57692 37.19013 36.98375 37.25068 
#>      L_1      L_3 
#> 36.49326 36.35529 
#>      H_1      H_2      H_3      L_1      L_2      H_1      H_2      H_3 
#> 18.48731 18.85805 19.58752 19.10595 19.55137 37.19013 36.98375 37.25068 
#>      L_1      L_2 
#> 36.49326 36.27569 
# }


## Using the SVD method with n_pcs set to 3.
imp_df4 <- impute_na(raw_df, method = "SVD", n_pcs = 3, seed = 3312)
#> change in estimate:  0.007642841 

## Using the minDet method with q set at 0.001.
imp_df5 <- impute_na(raw_df, method = "minDet", q = 0.001, seed = 3312)

## Impute a normalized data set using the kNN method
imp_df6 <- impute_na(ecoli_norm_df, method = "kNN")
#> Warning: Nothing to impute, because no NA are present (also after using makeNA)
```
