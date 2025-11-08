# Train machine learning models on training data

This function can be used to train models on protein intensity data
using different machine learning algorithms

## Usage

``` r
train_models(
  split_df,
  resample_method = "repeatedcv",
  resample_iterations = 10,
  num_repeats = 3,
  algorithm_list,
  seed = NULL,
  ...
)
```

## Arguments

- split_df:

  A `split_df` object from performing `split_data`.

- resample_method:

  The resampling method to use. Default is `"repeatedcv"` for repeated
  cross validation. See
  [`trainControl`](https://rdrr.io/pkg/caret/man/trainControl.html) for
  details on other available methods.

- resample_iterations:

  Number of resampling iterations. Default is `10`.

- num_repeats:

  The number of complete sets of folds to compute (For
  `resampling method = "repeatedcv"` only).

- algorithm_list:

  A list of classification or regression algorithms to use. A full list
  of machine learning algorithms available through the `caret` package
  can be found here:
  [http://topepo.github.io/caret/train-models-by-tag.html](http://topepo.github.io/caret/train-models-by-tag.md).
  See below for default options.

- seed:

  Numerical. Random number seed. Default is `NULL`

- ...:

  Additional arguments to be passed on to
  [`train`](https://rdrr.io/pkg/caret/man/train.html) function in the
  `caret` package.

## Value

A list of class `train` for each machine-learning algorithm. See
[`train`](https://rdrr.io/pkg/caret/man/train.html) for more information
on accessing different elements of this list.

## Details

- In the event that `algorithm_list` is not provided, a default list of
  four classification-based machine-learning algorithms will be used for
  building and training models. Default `algorithm_list`: "svmRadial",
  "rf", "glm", "xgbLinear, and "naive_bayes."

- Note: Models that fail to build are removed from the output.

- Make sure to fix the random number seed with `seed` for
  reproducibility

## References

Kuhn, Max. "Building predictive models in R using the caret package."
Journal of statistical software 28 (2008): 1-26.

## See also

- `pre_process`

- [`trainControl`](https://rdrr.io/pkg/caret/man/trainControl.html)

- [`train`](https://rdrr.io/pkg/caret/man/train.html)

## Author

Chathurani Ranathunge

## Examples

``` r
# \donttest{

## Create a model_df object
covid_model_df <- pre_process(covid_fit_df, covid_norm_df)
#> Total number of differentially expressed proteins (8) is less than n_top.
#> None of the proteins show high pair-wise correlation.
#> 
#> No highly correlated proteins to be removed.

## Split the data frame into training and test data sets
covid_split_df <- split_data(covid_model_df, seed = 8314)

## Fit models based on the default list of machine learning (ML) algorithms
covid_model_list1 <- train_models(split_df = covid_split_df, seed = 351)
#> 
#> Running svmRadial...
#> 
#> Running rf...
#> 
#> Running glm...
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> 
#> Running xgbLinear...
#> 
#> Running naive_bayes...
#> Done!

## Fit models using a user-specified list of ML algorithms.
covid_model_list2 <- train_models(
  covid_split_df,
  algorithm_list = c("svmRadial", "glmboost"),
  seed = 351
)
#> 
#> Running svmRadial...
#> 
#> Running glmboost...
#> glmboost failed.
#> Done!

## Change resampling method and resampling iterations.
covid_model_list3 <- train_models(
  covid_split_df,
  resample_method = "cv",
  resample_iterations = 50,
  seed = 351
)
#> Warning: `repeats` has no meaning for this resampling method.
#> 
#> Running svmRadial...
#> Warning: There were missing values in resampled performance measures.
#> 
#> Running rf...
#> Warning: There were missing values in resampled performance measures.
#> 
#> Running glm...
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: There were missing values in resampled performance measures.
#> 
#> Running xgbLinear...
#> Warning: There were missing values in resampled performance measures.
#> 
#> Running naive_bayes...
#> Warning: There were missing values in resampled performance measures.
#> Done!
# }
```
